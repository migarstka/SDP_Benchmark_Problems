# Code used to generate robust control problems.
# Save in Testproblems/RobustControl with an accurante
# solution from MOSEK


# workspace()
# include("./Helper.jl")

# using JuMP, Mosek, Base.Test, Helper, JLD
# Create number of Lyapunov stable controller and a number of Hinf performance controller
nnLyapunov = 40
nnHinf = 40

function LyapunovStabilityLMI(A,Bu)
  n = size(A,1)
  nu = size(Bu,2)
  AP  = kron(eye(n),A)
  PAt = kron(A,eye(n))
  BuY = kron(eye(n),Bu)
  YtBut = kron(Bu,eye(n))
  # turn YtBut*vec(Y') into YtBut* vec(Y) by permuting cols of YtBut
  YtBut = Helper.transposeVectorized(YtBut,n,nu) #provide dimension of Y'

  ba = spzeros(n^2,1)

  return [AP+PAt BuY+YtBut], ba
end


function HinfPerformanceLMI(A,Bu,Bw,Cz,Dzw,Dzu)
  n,nw = size(Bw)
  nz,nu = size(Dzu)
  # -----------------------------------
  # LMI): Hinf performance and stability
  # -----------------------------------

  # at first create the subblocks of the LMI left hand side
  # [A1     0       A2
  #  0      A3       0
  # A4      0       A5]
  # ---- A1 -----
  AP  = kron(eye(n),A)
  PAt = kron(A,eye(n))
  BuY = kron(eye(n),Bu)
  YtBut = kron(Bu,eye(n))
  YtBut = transposeVectorized(YtBut,n,nu) #provide dimension of Y'
  A1 = [AP+PAt BuY+YtBut zeros(n^2,1)] #will be multiplied by [vec(P),vec(Y),γ]^T
  # ---- A2 -----
  PCzt = kron(Cz,eye(n))
  YtDzut = kron(Dzu,eye(n))
  YtDzut = transposeVectorized(YtDzut,n,nu) #provide dimension of Y'
  A2 = [PCzt YtDzut zeros(n*nz,1)]
  # ---- A3 -----
  A3 = [zeros(nw^2,n^2+nu*n) -vec(eye(nw))]
  # ---- A4 -----
  CzP = kron(eye(n),Cz)
  DzuY = kron(eye(n),Dzu)
  A4 = [CzP DzuY zeros(nz*n,1)]
  # ---- A5 -----
  A5 = [zeros(nz^2,n^2+nu*n) -vec(eye(nz))]

  # Assemble big LMI block, s.t. Aa * [vec(P),vec(Y),γ]^T = vec(LMI)
  Aa = [A1[1:n,:];zeros(nw,size(A1,2));A4[1:nz,:]]
  for i = 1:n-1
    Aa = [Aa;A1[i*n+1:(i+1)*n,:];zeros(nw,size(A1,2));A4[i*nz+1:(i+1)*nz,:]]
  end
  for i = 0:nw-1
    Aa = [Aa;zeros(n,size(A1,2));A3[i*nw+1:(i+1)*nw,:];zeros(nz,size(A1,2))]
  end
  for i = 0:nz-1
    Aa = [Aa;A2[i*n+1:(i+1)*n,:];zeros(nw,size(A1,2));A5[i*nz+1:(i+1)*nz,:]]
  end

  Ba = [zeros(n,n) -Bw zeros(n,nz);
      -Bw' zeros(nw,nw) -Dzw';
      zeros(nz,n) -Dzw zeros(nz,nz)]

  ba = vec(Ba)
  # create right hand side matrix B, all constant entries of LMI
  return Aa,ba
end


rng = MersenneTwister(125)
dirPath = "../DataFiles/Julia/RobustControl/"

iii = 1
 while iii <=(nnLyapunov+nnHinf)

  # create random LTI system
  # choose number of states x(t)
  n = rand(rng,3:20)
  # choose number of outputs for z(t)
  nz = rand(rng,3:20)
  #choose number of inputs for u(t) and w(t)
  nu = rand(rng,3:20)
  nw = rand(rng,3:20)
  # create system matrices (create system in such a way that a stabilizing controller K exists)
  A = randn(rng,n,n)

  # choose Bw in such a way that (A,Bu) is controllable
  isControllable = false
  Bu = 0
  while !isControllable
    Bu = sprandn(rng,n,nu,0.8)
    # check controllability
    Con = Bu
    for kkk=1:n-1
      res = A^kkk*Bu
      Con = [Con res]
    end
    isControllable = (rank(full(Con)) == n)
  end
  Bw = eye(n,nw)
  Cz = randn(rng,nz,n)
  Dzw = randn(nz,nw)
  Dzu = eye(nz,nu)



  # LMI 1): Lyapunov stability: PA'+AP+Y'Bu'+BuY < 0 (NOTE STRICT INEQUALITY)
  if iii <=nnLyapunov
    Aa, ba = LyapunovStabilityLMI(A,Bu)
    Pa = spzeros(n^2+nu*n,n^2+nu*n)
    qa = spzeros(n^2+nu*n,1)
    nLMI = n
  else
    # LMI 2): Hinf performance:
    Aa,ba = HinfPerformanceLMI(A,Bu,Bw,Cz,Dzw,Dzu)
    nLMI = n+nw+nz
    Pa = spzeros(n^2+nu*n+1,n^2+nu*n+1)
    qa = spzeros(n^2+nu*n+1,1)
  end

  # # solve accurately once with mosek
  # model = Model(solver=MosekSolver())
  # @variable(model, P[1:n,1:n], SDP)
  # @variable(model, Y[1:nu,1:n])
  # @variable(model, g)
  # @objective(model, Min , g)
  # LMI = [A*P+P*A'+Bu*Y+Y'*Bu' Bw P*Cz'+Y'*Dzu';
  #         Bw' -g*eye(nw) Dzw';
  #         Cz*P+Dzu*Y Dzw -g*eye(nz)]
  # @SDconstraint(model, LMI <= 0)
  # status = JuMP.solve(model)
  # objTrue = getobjectivevalue(model)
  # P = getvalue(P)
  # Y = getvalue(Y)
  # g = getvalue(g)
  # K = Y/P
  # Ac = A+Bu*K
  # eig(Ac)

  # compute once gain in other format
  model2 = Model(solver=MosekSolver())
  @variable(model2, S2[1:nLMI,1:nLMI], SDP)
  @variable(model2, P2[1:n,1:n], SDP)
  @variable(model2, Y2[1:nu,1:n])
  if iii <= nnLyapunov
    @objective(model2, Min , 0)
    x2 = [vec(P2);vec(Y2)]
  else
    @variable(model2, g2)
    @objective(model2, Min , g2)
    x2 = [vec(P2);vec(Y2);g2]
  end
  s2 = vec(S2)
  @constraint(model2,Aa*x2+s2.==ba)
  status = JuMP.solve(model2)
  objTrue = getobjectivevalue(model2)
  P2 = getvalue(P2)
  Y2 = getvalue(Y2)
  K2 = Y2/P2
  Ac2 = A+Bu*K2
  # check eigenvalues of c.l. system A+Bu*F
  unstableK = maximum(map(x->real(x),eig(Ac2)[1])) >= 0.

  # discard all the weird cases with stalling solver or unstable eigenvalues
  if status != :Optimal || unstableK
    warn("Bad problem: Try again!")
    continue
  end

  # Some further tests for development
  # @test norm(Ac-Ac2,Inf) < 1e-5
  # @test abs(g-g2) < 1e-5
  # @test norm(vec(getvalue(LMI))-(Aa*[vec(P);vec(Y);g]-ba),Inf) < 1e-5
  # @test maximum(map(x->real(x),eig(Ac)[1])) <= 0.
  # @test maximum(map(x->real(x),eig(Ac2)[1])) <= 0.
  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end
  iii <= nnLyapunov ? pType = "LyapunovStability" : pType = "HinfPerformance"
  fn = "RobustControl"*nr*".jld"
  JLD.save(dirPath*fn,"n",size(Aa,2),"m",size(Aa,1),"A",Aa,"b",ba,"P",Pa,"q",qa,"objTrue",objTrue,"K",K2,"ProblemType",pType)
  println("$(iii)/$(nnLyapunov+nnHinf) completed!")
  iii+=1
end

