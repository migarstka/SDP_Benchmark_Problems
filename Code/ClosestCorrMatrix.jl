# Code used to generate closest correlation matrix problems. Save in Testproblems/ClosestCorr with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

using JuMP, Mosek, FileIO,SCS, Random, LinearAlgebra, SparseArrays, MathOptInterfaceMosek

# specify number of problems that start with matrix C as a perturbed correlation matrix and problems where C is a random matrix
#nnCorr = 0
#nnRand = 100
xMin = -1.
xMax = 1.
#nn = nnCorr+nnRand
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

# https://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-random-positive-semidefinite-correlation-matrices
# VINE METHOD to generate random correlation matrices
#distributed ~ det(S)^eta [or maybe det(S)^(eta-1), not sure]
function randomCorrMatrix(d, eta)
    β = eta + (d-1)/2;
    P = zeros(d,d);           #storing partial correlations
    S = eye(d);

    for k = 1:d-1
        β = β - 1/2;
        for i = k+1:d
            P[k,i] = rand(Beta(β,β)) # sampling from β
            P[k,i] = (P[k,i]-0.5)*2;     #linearly shifting to [-1, 1]
            p = P[k,i];
            for l = (k-1):-1:1 #converting partial correlation to raw correlation
                p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i]*P[l,k];
            end
            S[k,i] = p;
            S[i,k] = p;

        end
    end
    return S
end

rng = Random.MersenneTwister(12345);
dirPath = "../DataFiles/JLD2/SDP/ClosestCorr-Benchmark-v1/"
!ispath(dirPath) && mkdir(dirPath)

nRange = collect(25:25:800)
nn = length(nRange)
 for iii =1:1:nn
  # choose size of problem
  n = nRange[iii]
  eta = abs(randn(rng))


  C = xMin .+ rand(rng, n, n) .* (xMax - xMin)

  c = vec(C)
  n2 = n^2
  m = n+n2

  P = sparse(1.0I,n2,n2)
  q = -vec(C)
  r = 0.5*vec(C)'*vec(C)
  b = [ones(n);zeros(n2)]
  A = createDiagonalExtractor(n)
  Aa = [A; -sparse(1.0I,n2,n2)]

  # specify cone
  Kf = n
  Kl = 0
  Kq = []
  Ks = [n^2]

  # # # solve accurately once with mosek
  model = Model(with_optimizer(MosekOptimizer));
  @variable(model, X[1:n,1:n]);
  @variable(model, t);
  x = vec(X);
  @objective(model, Min, t);
  @constraint(model,  test,[t;x-c] in SecondOrderCone());
  @constraint(model, A*x.== b[1:n]);
  @SDconstraint(model, X >= 0);
  MOSEKtime = @elapsed status = JuMP.optimize!(model);

  # correct optimal objective value since slightly different problem is solved
  objTrue = JuMP.objective_value(model)
  MOSEKcost = 0.5*objTrue^2

  SCStime = 0.
  SCScost = 0.
  nr = "$(n)"
  if n < 100
    nr = "0$(n)"
  end

  fn = "ClosestCorr_N$(nr).jld2"
  problemType = "Closest Correlation Matrix"
  problemName = "ClosestCorr_N$(nr)"

  save(dirPath*fn,"n",n,"m",m,"A",Aa,"b",b,"P",P,"q",q,"r",r,"C",C,"problemType",problemType,"problemName",problemName,"Kf",Kf,"Ks",Ks,"MOSEKtime",MOSEKtime,"MOSEKcost",MOSEKcost,"SCStime",SCStime,"SCScost",SCScost)
  println("$(iii)/$(nn) completed!")
end


