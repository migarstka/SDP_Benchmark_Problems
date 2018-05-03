# Code used to generate SDP problems with a quadratic objective save in Testproblems/SDPQuad with an accurante
# solution from MOSEK
# The problem data is generated based on the Optimality conditions:
# Ax*+s*=b, s in S+, Px*+q+A'y*=0, y in S+

# Create different objective functions
# 1) P = rand pos def
# 2) P = I
# 3) P = 0

# A) Create problems like Ax+s=b, mat(S) in S+
# B) Create problems like <Ai,X> = bi, X in S+

# workspace()
# include("./Helper.jl")
# include("../../OSSDP/Code/src/Solver.jl")

# using JuMP, Mosek, JLD, HelperFunctions, OSSDP, Base.Test

rng = MersenneTwister(12345)
dirPath = "../DataFiles/Julia/SDPQuad/"
!ispath(dirPath) && mkdir(dirPath)
nn = 45

  # @testset "Random SDP Mosek vs. QOCS" begin

 for iii =1:1:nn
  # choose size of problem
  n = rand(rng,10:100)
  r = rand(rng,2:15)
  m = r^2
  A = spzeros(m,n)
  # construct A with several symmetric Ai that are vectorized and make the cols of A
  for kkk=1:n
    Ai = sprandn(rng,r,r,0.4)
    Ai = 0.5*(Ai+Ai')
    A[:,kkk] = vec(Ai)
  end
  xtrue = randn(rng,n,1)

  # create primal feasible point -> b
  Strue = generatePosDefMatrix(r,rng)
  strue = vec(Strue)
  b = A*xtrue + strue

  # create dual feasible point -> P,q
  # consider three different types of the quadratic cost
  if iii <= 15
     P = generatePosDefMatrix(n,rng)
     Phalf = P^(1/2)
  elseif iii <=30
     P = speye(n)
     Phalf = P
  else
   P = Phalf = spzeros(n,n)
  end

  Ytrue = generatePosDefMatrix(r,rng)
  ytrue = vec(Ytrue)
  q = (-P*xtrue -  A'*ytrue)[:]
  ra = 0.
  Kf = 0
  Kl = 0
  Kq = []
  Ks = [r^2]

  iii<=30 && (PhalfInvQ=Phalf\q)


  # K = Cone(0,0,[],[r^2])
  # setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
  # setON = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
  # setAdaptiveOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=true)
  # setAdaptiveON = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=true)
  # res1,nothing = OSSDP.solve(P,q,A,b,K,setOFF);
  # print("\n.")
  # res2,nothing = OSSDP.solve(P,q,A,b,K,setON);
  # print(".")
  # res3,nothing = OSSDP.solve(P,q,A,b,K,setAdaptiveOFF);
  # print(".")
  # res4,nothing = OSSDP.solve(P,q,A,b,K,setAdaptiveON);
  # print(".")

  # println(">>>QOCS: Problem $(iii), res1: cost: $(res1.cost), status: $(res1.status)")
  # println(">>>QOCS: Problem $(iii), res2: cost: $(res2.cost), status: $(res2.status)")
  # println(">>>QOCS: Problem $(iii), res3: cost: $(res3.cost), status: $(res3.status)")
  # println(">>>QOCS: Problem $(iii), res4: cost: $(res4.cost), status: $(res4.status)")

  # solve accurately once with mosek s
  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, S[1:r,1:r],SDP)
  # if P!= 0 add soc-constraint to model
  if P == spzeros(n,n)
    @objective(model, Min, q'*x)
  else
    # in case of quadratic objective the following problem is solved:
    # min x'Px + 2q'x
    # s.t. Ax+s=b, s in S+
    # instead solve equivalent problem
    # min t
    # s.t. ||Phalf*x+inv(Phalf)*q|| <= t
    # Ax+s=b, s in S+
    @variable(model, t)
    @objective(model, Min, t)
    @constraint(model,norm(Phalf*x+PhalfInvQ) <= t)
  end
  s = vec(S)
  @constraint(model,A*x+s.==b)
  status = JuMP.solve(model)
  objTrue = getobjectivevalue(model)

  (P!=spzeros(n,n)) && (objTrue = 0.5*( objTrue^2-q'*(P\q) ))
  solTrue = getvalue(x)

  # @test abs(objTrue-res1.cost) < 1e-2
  # @test norm(solTrue-res1.x,Inf) < 1e-2

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end




  fn = "SDPQuad"*nr*".jld"
  problemType = "Random SDP with quadratic Objective"
  problemName = "RandomSDP"*nr

  JLD.save(dirPath*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"r",ra,"objTrue",objTrue,"solTrue",solTrue,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks)
  println("$(iii)/$(nn) completed!")
end
  # end


