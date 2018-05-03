# Code used to generate SDP problems with a quadratic objective and a system of equality constraints save in Testproblems/SDPandEqConstraints with an accurante
# solution from MOSEK
# The problem data is generated based on the Optimality conditions:
# Ax*+s*=b, s in K, Px*+q+A'y*=0, y in K*

# steps:
# 1. Create random matrix A in R^mxn with 50% nonzeros
# 2. Half the elements in s belong to equality constraints, other half to one sdp cone, s=[s1;s2], with s1 in {0}-cone and s2 in psd-cone
# 3. Create random solution vector xtrue and strue with s1 = 0s, and mat(s2) in psd-cone
# 4. Determine right side of Ax+s=b
# 5. Create random pos def matrix P, and random dual variable ytrue=[y1;y2] with y1 in R+ and y2 in psd-cone
# 6. calculate q based on P,xtrue and ytrue

# workspace()
# include("./Helper.jl")
# include("../../OSSDP/Code/src/Solver.jl")

# using JuMP, Mosek, JLD, HelperFunctions, OSSDP, Base.Test

rng = MersenneTwister(12345)
dirPath = "../DataFiles/Julia/SDPandEqConstraints/"
!ispath(dirPath) && mkdir(dirPath)
nn = 30

   # @testset "Random SDP Mosek vs. QOCS" begin

 for iii =1:1:nn
  # choose size of problem
  n = rand(rng,10:100)
  r = rand(rng,2:10)
  # dimension of psd-cone (vector form)
  m2 = r^2
  # number of eq constraints
  m1 = m2
  m = m1+m2

  A = sprand(rng,m,n,0.5)
  xtrue = randn(rng,n,1)

  # create primal feasible point -> b
  s1 = zeros(m1,1)
  S2 = generatePosDefMatrix(r,rng)
  s2 = vec(S2)
  strue = [s1;s2]
  b = A*xtrue + strue

  # create dual feasible point -> P,q
   P = generatePosDefMatrix(n,rng)
   Phalf = P^(1/2)

  y1 = rand(rng,m1,1)
  Y2 = generatePosDefMatrix(r,rng)
  y2 = vec(Y2)
  ytrue = [y1;y2]
  q = (-P*xtrue -  A'*ytrue)[:]
  ra = 0.
  Kf = m1
  Kl = 0
  Kq = []
  Ks = [r^2]

  PhalfInvQ=Phalf\q


  # K = Cone(m1,0,[],[r^2])
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
  @variable(model, s1[1:m1])
  @variable(model, S2[1:r,1:r],SDP)

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
  s = [s1;vec(S2)]
  @constraint(model,A*x+s.==b)
  @constraint(model,s1.== 0.)
  status = JuMP.solve(model)
  objTrue = getobjectivevalue(model)

  objTrue = 0.5*( objTrue^2-q'*(P\q) )
  solTrue = getvalue(x)

 # @test abs(objTrue-res1.cost) < 1e-2
 # @test norm(solTrue-res1.x,Inf) < 1e-2

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "SDPandEqConstraints"*nr*".jld"
  problemType = "Random SDP with equality constraints"
  problemName = "SDPandEqConstraints"*nr

   JLD.save(dirPath*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"r",ra,"objTrue",objTrue,"solTrue",solTrue,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks)
  println("$(iii)/$(nn) completed!")
end
   # end


