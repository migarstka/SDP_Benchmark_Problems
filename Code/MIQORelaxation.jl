# Code used to generate Semidefinite Relaxation of Mixed Integer Quadratic Objective problems with an accurate
# solution from MOSEK

# Original problem has the following format:
# min_x   ||Ax-b||_2^2
# s.t.    x ∈ Z^n

workspace()
include("../../OSSDP/Code/src/Solver.jl")

using JLD, JuMP, Mosek, OSSDP


# this function creates a matrix A that slices out the diagonal entries Xii of a vectorized square matrix x=vec(X)
function createDiagonalExtractor(n)
  A = spzeros(n,n^2)
  A[1,1] = 1
  for iii=2:n-1
    col = (iii-1)*(n+1)
    A[iii,col+1] = 1
  end
  A[n,n^2] = 1
  return A
end

function createLMIconstMatrix(n)
  # assume the decision vector is like this [x,vec(X)]
  # M = spzeros((n+1)*(2*n),n+n^2)

  # # first add parts that belong to x and assign correct position in LMI
  # for iii=1:n
  #   E = spzeros(n+1,n+1)
  #   E[iii,end] = 1
  #   E[end,iii] = 1
  #   M[:,iii] = vec(E)
  # end

  # ind = n
  # # assign correct position of elements of X in LMI constraint
  # for jjj=1:n, iii=1:n
  #   E = spzeros(n+1,n+1)
  #   E[iii,jjj] = 1
  #   ind+=1
  #   M[:,ind] = vec(E)
  # end
  A1 = [zeros(n^2,n) kron(speye(n),speye(n))]

  # takes care of vec([X;x']) part
  Aa = [A1[1:n,:];1 zeros(1,n+n^2-1)]
  for i=1:n-1
    Aa = [Aa;A1[i*n+1:(i+1)*n,:];[zeros(1,i) 1 zeros(1,n+n^2-1-i)]]
  end
  # add vec([x;0]) part
  for i=0:n-1
    Aa = [Aa;zeros(1,i) 1 zeros(1,n+n^2-1-i)]
  end
  Aa = [Aa;spzeros(1,n+n^2)]
  return Aa
end

rng = MersenneTwister(12345)
dirPath =  "../DataFiles/Julia/MIQO/"
!ispath(dirPath) && mkdir(dirPath)

nn = 25


 for iii =1:1:nn
  # choose size of problem
  np = rand(rng,5:15)
  mp = 2*np
  A = randn(rng,mp,np)
  P = A'*A
  #xc: continous relaxation solution
  xc = rand(rng,np)
  q = -P*xc

  # scale problem to get fcts = -q'Pinv*q = -1
  # s = -q'*(P\q)
  # q = q./-s
  # P = P./-s
  # put problem into solver format
  # min   x'Ix - 1*t
  # s.t.  Ax + s1 = 0
  #       x = s4
  #       -t <= c'*x <=t
  # s1 in 0, s2,s3 in R+, s4 in SDPcone

  n = np+np^2
  m = (np+1)^2+np

  D = createDiagonalExtractor(np)
  LMI = createLMIconstMatrix(np)
  Aa = [speye(np,np) -D;-LMI]
  ba = [spzeros(np,1);spzeros((np+1)^2-1,1);1]
  Pa = spzeros(np+np^2,np+np^2)
  qa = [2*q;vec(P)]
  ra = 0

  # specify cone
  Kf = 0
  Kl = np
  Kq = []
  Ks = [(np+1)^2]

  # K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)
  # setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=5000,verbose=true,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
  # res1,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,setOFF);

  # model = Model()
  # @variable(model, 0 <= x[1:np] <= 1,Int)
  # @objective(model,Min, norm(A*x-b)^2)

  # status = JuMP.solve(model)


  # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, X[1:np,1:np])
  @variable(model, x[1:np])
  @objective(model,Min, trace(P*X) + 2*q'*x)

  @SDconstraint(model, [X x;x' 1] >= 0)
  for i=1:np
    @constraint(model, X[i,i] >= x[i])
  end
  status = JuMP.solve(model)

  objTrue = getobjectivevalue(model)
  solTrue = getvalue(x)


  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end
  problemType = "MIQO Relaxation"
  problemName = "MIQO"*nr

  fn = "MIQO"*nr*".jld"
  JLD.save(dirPath*fn,"n",n,"m",m,"A",Aa,"b",ba,"P",Pa,"q",qa,"r",ra,"objTrue",objTrue,"solTrue",solTrue,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks)
  println("$(iii)/$(nn) completed!")
end


