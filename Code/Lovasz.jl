# Code used to generate Lovasz Theta function problems. Save in Testproblems/Lovasz with an accurante
# solution from MOSEK

# Original problem has the following format:
# min_X  Tr(XJ)
# s.t.   Tr(X) = 1
#        X_ij = 1, if (i,j) ∈ E
#        X ⪴ 0

workspace()
using JLD, JuMP, Mosek



function createNonzeroExtractor(E)
  n = size(E,1)
  nz = nnz(sparse(E))
  S = spzeros(nz,n^2)
  e = vec(E)
  counter = 1
  for iii=1:length(e)
    if e[iii] == 1.
      S[counter,iii] = 1
      counter+=1
    end
  end
  return sparse(S)
end

# creates a random symmetric edge matrix with zeros on the diagonal and 1 as nonzero entries
function createRandomEdgeMatrix(rng,n,density)
  E  = full(Symmetric(sprand(rng,n,n,density)))
  map!(x->(x!=0 ? x=1. : x=0.),E),E
  for iii=1:n
    E[iii,iii] = 0
  end
  return sparse(E)
end

rng = MersenneTwister(12345)
dirPath = "./TestProblems/Lovasz/"
nn = 25


 for iii =1:1:nn
  # choose size of problem
  n = rand(rng,10:100)
  density = rand(rng,0.3:0.1:0.8)

  # create random edge matrix E for graph
  E = createRandomEdgeMatrix(rng,n,density)


  # Test: create E for pentagon graph
  # n = 5
  # E = zeros(5,5)
  # E[1,2] = E[2,3] = E[3,4] = E[4,5] = E[5,1] = 1
  # E[1,5] = E[2,1] = E[3,2] = E[4,3] = E[5,4] = 1
  # correctVal = sqrt(5)


  nz = nnz(sparse(E))
  J = sparse(ones(n,n))

  # put problem into solver format
  q = -vec(J)
  P = spzeros(n^2)

  A1 = vec(speye(n))'
  A2 = createNonzeroExtractor(E)
  A3 = -speye(n^2)

  Aa = [A1;A2;A3]
  ba = [1;zeros(nz+n^2)]


  # # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, X[1:n,1:n], SDP)
  @objective(model, Max ,trace(J*X))
  @constraint(model, trace(X) == 1)
  x = vec(X)
  @constraint(model, A2*x .== 0)
  status = JuMP.solve(model)

  objTrue = getobjectivevalue(model)
  solTrue = getvalue(x)
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


