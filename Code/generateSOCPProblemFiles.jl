# Test routine to compare scaling for a number of SOCP Lasso problems (partially badly scaled)
# Code used to generate QP Benchmark problems save in Testproblems/SOCP with an accurante
# solution from MOSEK
 # workspace()


 # using JuMP, Mosek, JLD

rng = MersenneTwister(12345)
dirPath = "./TestProblems/SOCP/"
nn = 25


for iii =1:1:nn
  # generate problem data
  n = 8
  m = 50*n
  F = rand(rng,m,n)
  vtrue = sprand(rng,n,1,0.1 )
  noise = 0.1*rand(rng,m,1)
  b = F*vtrue + noise
  μMax = norm(F'*b,Inf)
  μ = 0.1*μMax

  # define lasso problem as SOCP
    Aa = sparse([1 zeros(1,2*n+1) 1 zeros(1,m);
        -1 zeros(1,2*n) 1 zeros(1,m+1);
        zeros(m,1) -2*F zeros(m,n+2) eye(m,m);
        zeros(n,1) eye(n) -eye(n) zeros(n,m+2);
        zeros(n,1) -eye(n) -eye(n) zeros(n,m+2);
       zeros(1,2*n+1) -1 zeros(1,m+1);
       zeros(1,2*n+2) -1 zeros(1,m);
       zeros(m,2n+3) -eye(m)])

  ba = [1;1;-2*b;zeros(2*n+m+2)]
  q = [1;zeros(n);μ*ones(n,1);zeros(m+2,1)]
  P = spzeros(length(q),length(q))

  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, t[1:n])
  @variable(model, y)
  @variable(model, w[1:m+1])
  @objective(model, Min, y + μ*ones(n)'*t)
  @constraint(model, -t .<= x)
  @constraint(model, x .<= t)
  @constraint(model, w[1] == 1-y)
  @constraint(model, w[2:m+1] .== 2*(F*x-b))
  @constraint(model, norm(w) <= 1+y)
  status = JuMP.solve(model)
  println("Objective value: ", getobjectivevalue(model))
  println("x = ", getvalue(x))

  objTrue = getobjectivevalue(model)
  solTrue = getvalue(x)

  fn = "SOCP$(iii).jld"
  JLD.save(dirPath*fn,"n",n,"m",m,"A",sparse(Aa),"b",ba,"P",sparse(P),"q",q,"mu",μ,"objTrue",objTrue,"solTrue",solTrue)
  println("$(iii)/$(nn) completed!")
end


