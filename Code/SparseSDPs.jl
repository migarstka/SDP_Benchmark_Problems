include("./Helper.jl")
module SparseSDPs
using HelperFunctions
export generateArrowMultCones, generateBandedSDP

"""
generateSparseSymmetricMatrix(A::Array{Float64,2},amin::Int64,amax::Int64)
Returns a symmetric matrix based on the sparsity pattern of the input matrix A.
"""
function generateSparseSymmetricMatrix(rng,A::SparseMatrixCSC{Int64,Int64},aMin::Float64,aMax::Float64)
  m,n = size(A)
  m !=n && error("Input matrix must be square.")

  M = full(Symmetric(randn(rng,n,n)*(aMax-aMin)+aMin))
  for i = 2:n, j=1:i-1
    if A[i,j] == 0
     M[i,j] = 0.
     M[j,i] = 0.
   end
 end
 return sparse(M)
end

"""
generateArrowMultCones(rng,m,numCones,nBlk,BlkSize,ArrowWidth,NONZERO_P_FLAG)
Returns problem data P,q,A,b, Ks for a random problem with chordal sparsity in A and b.
"""
function generateArrowMultCones(rng,m,numCones,nBlk,BlkSize,ArrowWidth,NONZERO_P_FLAG)

  # sparsity pattern for each part of s is determined by respective rows in A and b
  Ks = zeros(Int64,numCones)

  # generate sparsity pattern based on user definition
  pattern = Array{SparseMatrixCSC{Int64,Int64}}(numCones)
  for iii=1:numCones
    dim = nBlk[iii]*BlkSize[iii]+ArrowWidth[iii]
    pattern[iii] = spzeros(dim,dim)
    p = pattern[iii]
    for k=1:nBlk[iii]
      p[(k-1)*BlkSize[iii]+1:BlkSize[iii]*k, (k-1)*BlkSize[iii]+1:BlkSize[iii]*k] = 1
    end
    # create the arrow head
    p[nBlk[iii]*BlkSize[iii]+1:dim,:] = 1
    p[:,nBlk[iii]*BlkSize[iii]+1:dim] = 1
    # save cone size
    Ks[iii] = dim

  end

  numRows = sum(map(x->x^2,Ks))
  A = spzeros(numRows,m)
  # generate random problem data based on pattern
  for iii=1:m
    Ai = Array{Float64}(0,1)
    for k= 1:numCones
      M = generateSparseSymmetricMatrix(rng,pattern[k],-25.,25.)
      Ai = [Ai;vec(M)]
    end
    A[:,iii] = Ai
  end

  # create feasible primal point
  xtrue = rand(rng,m)

  strue = spzeros(size(A,1))
  b = 1
  for k=1:numCones
    Temp = generateSparseSymmetricMatrix(rng,pattern[k],-10.,10.)
    dim = size(pattern[k],1)
    Temp = Temp + (-minimum(eigs(Temp)[1])+1)*speye(dim)
    e = Int(b + dim^2) - 1
    strue[b:e] = vec(Temp)
    b = e+1
  end
  ba = sparse(A*xtrue + strue)

  # generate feasible dual point ytrue
  if NONZERO_P_FLAG
    P = HelperFunctions.generatePosDefMatrix(m,rng,0.1,5)
  else
    P = spzeros(m,m)
  end

  ytrue = Array{Float64}(0)
  for k=1:numCones
    Ytrue = generatePosDefMatrix(Ks[k],rng,0.1,5)
    ytrue = [ytrue;vec(Ytrue)]
  end
  q = (-P*xtrue -  A'*ytrue)[:]

  map!(x->x^2,Ks,Ks)
  return P,q,A,full(ba),Ks
end


function randomBandedSymMatrix(rng,n::Int64,bw::Int64,aMin::Float64,aMax::Float64)
  A = spzeros(n,n)

  for i=0:bw
    r = rand(rng,n-i)*(aMax-aMin)+aMin
    A = A+diagm(r,i)
    i != 0 && (A = A+diagm(r,-i))
  end
  return A
end


function generateBandedSDP(rng::MersenneTwister,m::Int64,K,bandwidth::Array{Int64},aMin::Float64,aMax::Float64,density::Float64)

  for iii=1:length(bandwidth)
    bw = bandwidth[iii]
    bw >= Int(sqrt(K.s[iii])) && error("Your bandwidth is larger than you problem matrix.")
  end


  numRows = K.f + K.l + sum(K.q) + sum(K.s)
  A = spzeros(numRows,m)
  for i = 1:m
    Ai = sprand(rng,K.f,1,density)
    Ai = [Ai;sprand(rng,K.l,1,density)]
    Ai = [Ai;sprand(rng,sum(K.q),1,density)]
    # loop over cones and create banded structure
    for k=1:length(K.s)
      dim = Int(sqrt(K.s[k]))
      Temp = randomBandedSymMatrix(rng,dim,bandwidth[k],aMin,aMax)
      Ai = [Ai;vec(Temp)]
    end
    A[:,i] = Ai
  end

  # create strictly feasible primal point
  xtrue = rand(rng,m)

  strue = spzeros(K.f)
  strue = [strue;rand(rng,K.l)*10]
  for i=1:length(K.q)
    Temp = rand(rng,K.q[i]-1)
    strue = [strue;norm(Temp,2)+1;Temp]
  end

  for i=1:length(K.s)
    dim = Int(sqrt(K.s[i]))
    Temp = randomBandedSymMatrix(rng,dim,bandwidth[i],aMin,aMax)
    Temp = Temp + (-minimum(eigs(Temp)[1])+1)*speye(dim)
    strue = [strue;vec(Temp)]
  end
  b = A*xtrue + strue



  # create strictly feasible dual point
  P = generatePosDefMatrix(m,rng)

  ytrue = rand(rng,K.f+K.l)*10
  for i=1:length(K.q)
    Temp = rand(rng,K.q[i]-1)
    ytrue = [ytrue;norm(Temp,2)+1;Temp]
  end

  for i=1:length(K.s)
    dim = Int(sqrt(K.s[i]))
    Ytrue = generatePosDefMatrix(dim,rng)
    ytrue = [ytrue;vec(Ytrue)]
  end
  q = (-P*xtrue -  A'*ytrue)[:]

  return P,q,A,b,K
end


end #END Module