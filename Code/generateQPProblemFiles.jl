# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)
# workspace()
# Code used to generate QP Benchmark problems save in Testproblems/QP with an accurante
# solution from MOSEK

 using JuMP, Mosek, JLD

rng = MersenneTwister(12345)
dirPath = "./TestProblems/QP/"
nn = 25

for iii =1:1:nn
  # generate problem data
  n = 8
  m = 50*n
  A = sprandn(rng,m,n,0.5)
  vtrue = 1/n*sprandn(rng,n,0.5)
  noise = 1/4*randn(rng,m,1)
  b = A*vtrue + noise
  λ = 0.2*norm(A'*b,Inf)

  pDims = [m;n;nnz(A)]

  # define lasso problem as QP
  Aa = [-A zeros(m,n) eye(m,m);
         eye(n,n) -eye(n,n) zeros(n,m);
         -eye(n,n) -eye(n,n) zeros(n,m)]

  ba = [-b;zeros(2*n)]
  P = 2*diagm([zeros(2*n);ones(m)])
  q = [zeros(n);λ*ones(n);zeros(m)]

  # modify problem for OSQP (using inequality constraints)
  l_OSQP = full([-b;-Inf*ones(n);zeros(n)][:])
  u_OSQP = full([-b;zeros(n);Inf*ones(n)][:])
  A_OSQP = [-A zeros(m,n) eye(m,m);
       eye(n,n) -eye(n,n) zeros(n,m);
       eye(n,n) eye(n,n) zeros(n,m)]

  # solve with MOSEK for accurate solution
  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, t[1:n])
  @variable(model, y[1:m])
  @objective(model, Min, y'y + λ*ones(n)'*t)
  @constraint(model, -t .<= x)
  @constraint(model, x .<= t)
  @constraint(model, y .== A*x-b)
  status = JuMP.solve(model)

  objTrue = getobjectivevalue(model)
  solTrue = getvalue(x)

  fn = "QP$(iii).jld"
  JLD.save(dirPath*fn,"n",n,"m",m,"A",sparse(Aa),"b",ba,"P",sparse(P),"q",q,"lambda",λ,"objTrue",objTrue,"solTrue",solTrue,"A_OSQP",A_OSQP,"l_OSQP",l_OSQP,"u_OSQP",u_OSQP)
  println("$(iii)/$(nn) completed!")
end


