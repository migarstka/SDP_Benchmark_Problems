# Code used to generate smallest sphere around multiple ellipsoids problems. Save in Testproblems/SmallestSphere with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_{t,xc,γ,τ1,...,τk)   t
# s.t.    [I -xc; xc' t+γ] ⪳ τi *  [Ai bi; bi' ci], i+1,...,k
#         [I, xc; xc' t+γ] ⪴ 0
#         τi ≥ 0, i+1,...,k

workspace()
#include("../../src/Solver.jl")
using JuMP, Mosek, JLD
#using PyPlot, PyCall
#@pyimport matplotlib.patches as patch

rng = MersenneTwister(12345)
dirPath = "./TestProblems/SmallestSphere/"
nn = 25


 for iii =1:1:nn
  # set dimension of ellipsoids
  n = 3
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
    Peigs = rand(rng,0.25:0.01:8,3)
    α = rand(rng,-pi/2:0.01:pi/2)
    β = rand(rng,-pi/2:0.01:pi/2)
    ϕ = rand(rng,-pi/2:0.01:pi/2)
    P = diagm(Peigs)
    Rx = [1 0 0;0 cos(α) -sin(α); sin(α) cos(α)]
    Ry = [cos(β) 0 -sin(β);0 1 0; -sin(β) 0 cos(β)]
    Rz = [cos(ϕ) -sin(ϕ) 0; sin(ϕ) cos(ϕ) 0; 0 0 1]
    R = Rx*Ry*Rz
    P = R*P*R'
    xc = rand(rng,-10:10,3,1)

    # determine corresponding matrices for sublevel representation f(x)=x'Ax+2x'b+c
    Ainv = P\eye(3)
    A = Ainv*Ainv
    b = -A*xc
    c = (xc'A*xc-1)[1]
    push!(Ae,A)
    push!(be,b)
    push!(ce,c)
    push!(Pe,P)
    push!(xce,xc)
  end

  # put problem into solver format -> P,q,Aa,ba,K
  # Define order of decision variables x=[t;xc;γ;τ1;...;τk]
  dimLMI = n+1

  # 1) constraint 1: τi ≥ 0, i+1,...,k
  A1 = [spzeros(ne,2+n) -speye(ne)]
  Aa = A1
  ba = zeros(ne,1)

  # 2) constraint 2: [I, xc; xc' t+γ] ⪴ 0
  # NOTE: HERE I ASSUME DIM(xc) = 2
  Mt = zeros(dimLMI,dimLMI)
  Mt[dimLMI,dimLMI] = 1
  Mγ = Mt
  Mxc1 = zeros(dimLMI,dimLMI)
  Mxc1[1,3] = 1
  Mxc1[3,1] = 1
  Mxc2 = zeros(dimLMI,dimLMI)
  Mxc2[2,3] = 1
  Mxc2[3,2] = 1
  M0 = zeros(dimLMI,dimLMI)
  M0[1:n,1:n] = eye(n)
  A2 = [vec(Mt) vec(Mxc1) vec(Mxc2) vec(Mγ) zeros(dimLMI^2,ne)]

  # append LMI constraint to Augmented Data Matrix Aa and ba
  Aa = [Aa;-A2]
  ba = [ba;vec(M0)]

  # 3) constraint 3: [I -xc; xc' t+γ]  ⪳ τi *  [Ai bi; bi' ci], i+1,...,k
  Mt = zeros(dimLMI,dimLMI)
  Mγ = zeros(dimLMI,dimLMI)
  Mγ[dimLMI,dimLMI] = 1
  Mxc1 = zeros(dimLMI,dimLMI)
  Mxc1[1,3] = -1
  Mxc1[3,1] = -1
  Mxc2 = zeros(dimLMI,dimLMI)
  Mxc2[2,3] = -1
  Mxc2[3,2] = -1
  M0 = zeros(dimLMI,dimLMI)
  M0[1:n,1:n] = eye(n)
  for kkk=1:ne
    Mτ = zeros(dimLMI,dimLMI)
    Mτ[1:n,1:n] = -Ae[kkk]
    Mτ[1:n,n+1] = -be[kkk]
    Mτ[n+1,1:n] = -be[kkk]'
    Mτ[dimLMI,dimLMI] = -ce[kkk]

    # append LMI constraint to Augmented Data Matrix Aa and ba
    A3 = [vec(Mt) vec(Mxc1) vec(Mxc2) vec(Mγ) zeros(dimLMI^2,kkk-1) vec(Mτ) zeros(dimLMI^2,ne-kkk)]
    Aa = [Aa;A3]
    ba = [ba;-vec(M0)]
  end

  Pa = spzeros(size(Aa,2),size(Aa,2))
  qa = [1;zeros(size(Aa,2)-1)]

  # K = OSSDPTypes.Cone(0,ne,[],dimLMI^2*ones(ne+1))
  # setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=5000,verbose=true,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
  # res1,nothing = OSSDP.solve(Pa,qa,Aa,ba,K,setOFF)
  # xc = res1.x[2:3]
  # γ = res1.x[4]
  # r = sqrt(xc'*xc-γ)


  # # plot ellipsoids

  fig = PyPlot.figure(1,facecolor="white",figsize=(12,5))
  ax = fig[:add_subplot](1,1,1)
  ax[:set_aspect]("equal")

  for i = 1:ne
    # convert before plotting
      xs = Float64[]
      ys = Float64[]
      for angle in linspace(0, 2*pi, 100)
          u = [cos(angle),sin(angle)]
          x = Pe[i] * u
          x += xce[i]
          push!(xs, x[1])
          push!(ys, x[2])
      end
      PyPlot.plot(xs, ys, "black", linewidth=1.0)
  end
  # plot solution sphere
   c = patch.Circle([xc[1],xc[2]],r,ec="red",fill=false, linewidth=1.,zorder=0)
   ax[:add_artist](c)
  PyPlot.axis([-20 20 -20 20]')

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

  fn = "SmallestCircle"*nr*".jld"
  JLD.save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",Pa,"q",qa,"objTrue",objTrue,"r",r,"xc",xc)
  println("$(iii)/$(nn) completed!")
end


