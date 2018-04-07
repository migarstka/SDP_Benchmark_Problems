# Code used to generate polynomial sum-of-squares problems. Save in Testproblems/PolynomialSOS with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_Q   0
# s.t.    p(x) = x'Qx,
#         Q ⪴ 0

workspace()
#include("../../src/Solver.jl")
using JuMP, Mosek, JLD

rng = MersenneTwister(12345)
dirPath = "./TestProblems/PolySOS/"
nn = 25


 for iii =1:1:nn
  # set dimension of ellipsoids
  n = 2
  # choose number of elipsoids
  ne = rand(rng,2:10)
  # create random ellipsoid data
  Ae = []
  be = []
  ce = []
  Pe = []
  xce = []
  for kkk=1:ne
    # create mapping matrix P that distorts unit circle by scaling principal axis and rotating
    # define center of ellisoid xc
    Peigs = rand(rng,0.25:0.01:8,2)
    theta = rand(rng,-pi:0.01:pi)
    P = diagm(Peigs)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    P = R*P*R'
    xc = rand(rng,-10:10,2,1)

    # determine corresponding matrices for sublevel representation f(x)=x'Ax+2x'b+c
    Ainv = P\eye(2)
    A = Ainv*Ainv
    b = -A*xc
    c = (xc'A*xc-1)[1]
    push!(Ae,A)
    push!(be,b)
    push!(ce,c)
    push!(Pe,P)
    push!(xce,xc)
  end

  #solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, xc[1:2])
  @variable(model, t)
  @variable(model, g)
  @variable(model, tau1 >= 0)
  @variable(model, tau2 >= 0)

  @objective(model, Min, t)
  @SDconstraint(model,[eye(2) xc;xc' t+g] >= 0)
  @SDconstraint(model,[eye(2) -xc;-xc' g] <= tau1*[Ae[1] be[1]; be[1]' ce[1]])
  @SDconstraint(model,[eye(2) -xc;-xc' g] <= tau2*[Ae[2] be[2]; be[2]' ce[2]])
  status = JuMP.solve(model)

  # correct optimal objective value since slightly different problem is solved
  objTrue = getobjectivevalue(model)
  xc = getvalue(xc)
  g = getvalue(g)
  r = sqrt(xc'*xc-g)

  # println("Objective value: ",objTrue)
  # println("x = ", getvalue(x))

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "PolySOS"*nr*".jld"
  JLD.save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",Pa,"q",qa,"objTrue",objTrue,"r",r,"xc",xc)
  println("$(iii)/$(nn) completed!")
end


