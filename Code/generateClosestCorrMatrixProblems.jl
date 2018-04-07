# Code used to generate closest correlation matrix problems. Save in Testproblems/ClosestCorr with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

workspace()
using JLD, JuMP, Mosek


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

rng = MersenneTwister(12345)
dirPath = "./TestProblems/ClosestCorr/"
nn = 25


 for iii =1:1:nn
  # choose size of problem
  n = rand(rng,6:40)
  eta = abs(randn(rng))
  # create random correlation matrix and perturb it
  C = randomCorrMatrix(n,eta)
  E  = randn(n,n)*1e-1
  C = C + (E + E')/2
  c = vec(C)

  isposdef(C) && warn("The perturbed correlation matrix is still pos def.")

  # put problem into solver format
  # min   x'Ix - 1*t
  # s.t.  Ax + s1 = 0
  #       x = s4
  #       -t <= c'*x <=t
  # s1 in 0, s2,s3 in R+, s4 in SDPcone

  n2 = n^2
  m = n+n2

  P = speye(n2)
  q = -vec(C)
  b = [ones(n);zeros(n2)]
  A = createDiagonalExtractor(n)
  Aa = [A; -speye(n2)]

  # # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, X[1:n,1:n], SDP)
  @variable(model, t)
  x = vec(X)
  @objective(model, Min, t)
  @constraint(model, norm(x-c) <= t)
  @constraint(model, A*x.== b[1:n])
  status = JuMP.solve(model)

  # correct optimal objective value since slightly different problem is solved
  objTrue = getobjectivevalue(model)
  objCorrected = 0.5*objTrue^2
  solTrue = getvalue(x)
  # println("Objective value: ",objTrue)
  # println("Corrected objective value: ", objCorrected)
  # println("x = ", getvalue(x))

  nr = "$(iii)"
  if iii < 10
    nr = "0$(iii)"
  end

  fn = "ClosestCorr"*nr*".jld"
  JLD.save(dirPath*fn,"n",n,"m",m,"A",Aa,"b",b,"P",P,"q",q,"C",C,"objTrue",objCorrected,"solTrue",solTrue)
  println("$(iii)/$(nn) completed!")
end


