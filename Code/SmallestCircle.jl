# Code used to generate smallest circle around multiple ellipsoids problems. Save in Testproblems/SmallestCircle with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_{t,xc,γ,τ1,...,τk)   t
# s.t.    [I -xc; xc' t+γ] ⪳ τi *  [Ai bi; bi' ci], i+1,...,k
#         [I, xc; xc' t+γ] ⪴ 0
#         τi ≥ 0, i+1,...,k

workspace()
include("../../OSSDP/Code/src/Solver.jl")
using JuMP, Mosek, JLD


rng = MersenneTwister(12345)
dirPath = "../DataFiles/Julia/SmallestCircle/"
!ispath(dirPath) && mkdir(dirPath)
nn = 25

 for iii =1:1:nn
  # set dimension of ellipsoids
  n = 2
  # choose number of elipsoids
  ne = rand(rng,3:15)
  # create random ellipsoid data
  Ae = []
  be = []
  ce = []
  Pe = []
  xce = []
  for kkk=1:ne
    # create mapping matrix P that distorts unit circle by scaling principal axis and rotating
    # define center of ellisoid xc
    Peigs = rand(rng,2)*(8-0.25)+0.25
    theta = rand(rng)*2*pi-pi
    P = diagm(Peigs)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    P = R*P*R'
    xc = rand(rng,2,1)*(10.-(-10))-10

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

  # put problem into solver format -> Pa,qa,Aa,ba,K
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

  # append LMI constraint to augmented Data Matrix Aa and ba
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
  ra = 0.

  # define cone dimensions
  Kf = 0
  Kl = ne
  Kq = []
  Ks = dimLMI^2*ones(ne+1)



  # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, xc[1:2])
  @variable(model, t)
  @variable(model, g)
  @variable(model, tau[1:ne] >= 0)
  @objective(model, Min, t)
  @SDconstraint(model,[eye(2) xc;xc' t+g] >= 0)
  for jjj=1:ne
    @SDconstraint(model,[eye(2) -xc;-xc' g] <= tau[jjj]*[Ae[jjj] be[jjj]; be[jjj]' ce[jjj]])
  end
  status = JuMP.solve(model)

  @test status == :Optimal
  objTrue = getobjectivevalue(model)
  xc = getvalue(xc)
  g = getvalue(g)
  radius = sqrt(xc'*xc-g)

  # # try mosek with qocs solver format
  # model2 = Model(solver=MosekSolver())
  # @variable(model2, S0[1:ne] >= 0)
  # @variable(model2, S1[1:dimLMI,1:dimLMI],SDP)
  # @variable(model2, S2[1:dimLMI,1:dimLMI,1:ne])
  # @variable(model2, xc2[1:2])
  # @variable(model2, t2)
  # @variable(model2, g2)
  # @variable(model2, tau2[1:ne] >= 0)

  # x = [t2;xc2;g2;tau2]ss
  # s = [S0;vec(S1)]
  # for jjj=1:ne
  #   s = [s;vec(S2[:,:,jjj])]
  #   @SDconstraint(model2,S2[:,:,jjj] >= 0)
  # end
  # @objective(model2, Min, qa'*x)
  # @constraint(model2, Aa*x+s.==ba)
  # status2 = JuMP.solve(model2)

  # objTrue2 = getobjectivevalue(model2)
  # xc2 = getvalue(xc2)
  # g2 = getvalue(g2)
  # radius2 = sqrt(xc2'*xc2-g2)
  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end
  problemType = "Smallest Circle around ellipsoids"
  problemName = "SmallestCircle"*nr
  fn = "SmallestCircle"*nr*".jld"
  JLD.save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",Pa,"q",qa,"r",ra,"objTrue",objTrue,"radius",radius,"xc",xc,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks)
  println("$(iii)/$(nn) completed!")
end

end
