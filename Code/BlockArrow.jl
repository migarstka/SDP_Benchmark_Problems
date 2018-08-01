# Code used to generate random sdp problems with block arrow pattern
workspace()
include("./SparseSDPs.jl")
 using JLD, SparseSDPs, JuMP, Mosek


rng = MersenneTwister(13123)
dirPath = "../DataFiles/DecomposableProblems/BlockArrow/"
NONZERO_P_FLAG = false

mRange = collect(5:5:75)
lRange = collect(3:3:45)
dRange = collect(2:1:16)


function solveWithMOSEK(P,q,A,b,Ks)
  m,n = size(A)
  cDim = Int(sqrt(m))
  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, S[1:cDim,1:cDim],SDP)
  s = vec(S)
  @objective(model, Min,q'*x)
  @constraint(model, A*x.+s .== b)
  tic()
  status = JuMP.solve(model)
  solveTime = toq()
  return getobjectivevalue(model),solveTime
end


# --------------------------------------
# Variable number of constraints m
# --------------------------------------
println(">>Start creating problems with variable number of constraints m:")
 for iii =1:1:length(mRange)
  nn = length(mRange)
  numCones = 1
  numBlocks = [10]
  BlkSize = [4]
  m = mRange[iii]
  ArrowWidth = [4]

  P,q,A,b,Ks = generateArrowMultCones(rng, m,numCones,numBlocks,BlkSize,ArrowWidth,NONZERO_P_FLAG);
  Kf  = 0
  Kl = 0
  Kq = []
  r = 0.
  m,n = size(A)
  objTrue,solveTime = solveWithMOSEK(P,q,A,b,Ks)
  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "BlkArrow_varM"*nr*".jld"
  problemType = "BlkArrow_varM"
  problemName = "BlkArrow_varM"*nr
  extraDir = "varM/"
  !ispath(dirPath*extraDir) && mkdir(dirPath*extraDir)

  JLD.save(dirPath*extraDir*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"r",r,"l",numBlocks[1],"d",BlkSize[1],"mconstr",mRange[iii],"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks,"objTrue",objTrue,"solveTime",solveTime)
  println("$(iii)/$(nn) completed!")
end

# --------------------------------------
# Variable number of number of blocks l
# --------------------------------------
println(">>Start creating problems with variable number of blocks l:")

# variable constraint number problems
 for iii =1:1:length(lRange)
  nn = length(lRange)

  numCones = 1
  numBlocks = lRange[iii]
  BlkSize = [4]
  m = 20
  ArrowWidth = [4]

  P,q,A,b,Ks = generateArrowMultCones(rng, m,numCones,numBlocks,BlkSize,ArrowWidth,NONZERO_P_FLAG);
  Kf  = 0
  Kl = 0
  Kq = []
  r = 0.
  m,n = size(A)
  objTrue,solveTime = solveWithMOSEK(P,q,A,b,Ks)

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "BlkArrow_varL"*nr*".jld"
  problemType = "BlkArrow_varL"
  problemName = "BlkArrow_varL"*nr
  extraDir = "varL/"
  !ispath(dirPath*extraDir) && mkdir(dirPath*extraDir)

  JLD.save(dirPath*extraDir*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"r",r,"l",lRange[iii],"d",BlkSize[1],"mconstr",20,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks,"objTrue",objTrue,"solveTime",solveTime)
  println("$(iii)/$(nn) completed!")
end

# --------------------------------------
# Variable Block size d
# --------------------------------------
println(">>Start creating problems with variable block size d:")

# variable constraint number problems
 for iii =1:1:length(dRange)
  nn = length(dRange)
  numCones = 1
  numBlocks = [10]
  BlkSize = [dRange[iii]]
  m = 10
  ArrowWidth = [4]

  P,q,A,b,Ks = generateArrowMultCones(rng, m,numCones,numBlocks,BlkSize,ArrowWidth,NONZERO_P_FLAG);
  Kf  = 0
  Kl = 0
  Kq = []
  r = 0.
  m,n = size(A)
  objTrue,solveTime = solveWithMOSEK(P,q,A,b,Ks)

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "BlkArrow_varD"*nr*".jld"
  problemType = "BlkArrow_varD"
  problemName = "BlkArrow_varD"*nr
  extraDir = "varD/"
  !ispath(dirPath*extraDir) && mkdir(dirPath*extraDir)

  JLD.save(dirPath*extraDir*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"r",r,"l",numBlocks[1],"d",dRange[iii],"mconstr",10,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Kl",Kl,"Kq",Kq,"Ks",Ks,"objTrue",objTrue,"solveTime",solveTime)
  println("$(iii)/$(nn) completed!")
end


