# Test routine to compare scaling for a number of QP Lasso problems (partially badly scaled)

using Test, LinearAlgebra, SparseArrays, Random, Statistics, FileIO

rng = MersenneTwister(12345)

dirPath = "../DataFiles/JLD2/QP/QP-Lasso/"
!ispath(dirPath) && mkdir(dirPath)

nrange = collect(20:10:220)
nn = length(nrange)
for iii =1:1:length(nrange)
  # generate problem data
  n = nrange[iii]
  m = 50*n
  A = sprandn(rng,m,n,0.5)
  vtrue = 1/n*sprandn(rng,n,0.5)
  noise = 1/4*randn(rng,m,1)
  b = A*vtrue + noise
  λ = 0.2*norm(A'*b,Inf)


  # define lasso problem as QP
  Aa = [-A zeros(m,n) Matrix(1.0I,m,m);
         Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m);
         -Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m)]

  ba = [-b;zeros(2*n)]
  P = sparse(2*Diagonal([zeros(2*n);ones(m)])) # times two to cancel the 1/2 in the cost function
  q = [zeros(n);λ*ones(n);zeros(m)]
  r = 0.
  # define cone membership
  Kf = m
  Kl = 2*n
  Kq = []
  Ks = []


  # modify problem for OSQP (using inequality constraints)
  l = Vector([-b;-Inf*ones(n);zeros(n)][:])
  u = Vector([-b;zeros(n);Inf*ones(n)][:])
  A_OSQP = [-A zeros(m,n) Matrix(1.0I,m,m);
       Matrix(1.0I,n,n) -Matrix(1.0I,n,n) zeros(n,m);
       Matrix(1.0I,n,n) Matrix(1.0I,n,n) zeros(n,m)]

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "QL-Lasso"*nr*".jld2"
  problemType = "QP-Lasso"
  problemName = "QP-Lasso"*nr

  save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",P,"q",q,"r",r,"A_OSQP",A_OSQP,"l_OSQP",l,"u_OSQP",u,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks)
  println("$(iii)/$(nn) completed!")
end


