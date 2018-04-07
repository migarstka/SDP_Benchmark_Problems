# Code used to generate robust control problems. Save in Testproblems/RobustControl with an accurante
# solution from MOSEK


workspace()
using JLD, JuMP, Mosek




rng = MersenneTwister(12345)
dirPath = "./TestProblems/RobustControl/"
nn = 25


 for iii =1:1:nn
  # choose number of states x(t)
  n = rand(rng,3:5)
  # choose number of outputs for z(t)
  nz = n
  #choose number of inputs for u(t) and w(t)
  nu = rand(rng,1:3)
  nw = rand(rng,1:3)
  # create system matrices (create system in such a way that a stabilizing controller K exists)
  A = sprandn(n,n,0.8)

  # choose Bw in such a way that (A,Bu) is controllable
  isControllable = false
  Bu = 0
  while !isControllable
    Bu = sprandn(n,nu,0.8)
    # check controllability
    Con = Bu
    for kkk=1:n-1
      res = A^kkk*Bu
      Con = [Con res]
    end
    isControllable = (rank(full(Con)) == n)
  end
  Bw = zeros(n,nw)
  Cz = eye(nz)
  Dzw = eye(nz,nw)
  Dzu = zeros(nz,nu)

  # LMI 1): Lyapunov stability: PA'+AP+Y'Bu'+BuY < 0 (NOTE STRICT INEQUALITY)


  # # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, P[1:n,1:n], SDP)
  @variable(model, Y[1:nu,1:n])
  @variable(model, g)
  @objective(model, Min , g)
  LMI = [A*P+P*A'+Bu*Y+Y'*Bu' Bw P*Cz'+Y'*Dzu';
          Bw' -g*eye(nw) Dzw';
          Cz*P+Dzu*Y Dzw -g*eye(nz)]
  @SDconstraint(model, LMI <= 0)
  status = JuMP.solve(model)

  objTrue = getobjectivevalue(model)
  P = getvalue(P)
  Y = getvalue(Y)
  F = Y/P

  # compute closed loop system matrix Ac
  Ac = A+Bu*F
  eig(Ac)
  # compute Hinf gain


  # println("Objective value: ",objTrue)
  # println("Corrected objective value: ", objCorrected)
  # println("x = ", getvalue(x))

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "Lovasz"*nr*".jld"
  JLD.save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",P,"q",q,"E",E,"objTrue",objTrue,"solTrue",solTrue)
  println("$(iii)/$(nn) completed!")
end


