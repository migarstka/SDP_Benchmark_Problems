# Code used to generate SDP problems with a quadratic function save in Testproblems/SDPQuad with an accurante
# solution from MOSEK
# The problem data is generated based on the Optimality conditions:
# Ax*+s*=b, s* in S+, Px*+q+A'y*=0, y in S+

# Create different objective functions
# 1) P = rand pos def
# 2) P = I
# 3) P = 0

# A) Create problems like Ax+s=b, mat(S) in S+
# B) Create problems like <Ai,X> = bi, X in S+

workspace()
include("../../src/Helper.jl")
using JuMP, Mosek, JLD, Helper

rng = MersenneTwister(12345)
dirPath = "./TestProblems/SDPQuad/"
nn = 25


# for iii =1:1:nn
  # choose size of problem
  n = rand(rng,10:100)
  r = rand(rng,2:10)
  m = r^2
  A = zeros(m,n)
  # construct A with several symmetric Ai that are vectorized and make the cols of A
  for kkk=1:n
    Ai = sprandn(rng,r,r,0.4)
    Ai = 0.5*(Ai+Ai')
    A[:,kkk] = vec(Ai)
  end
  xtrue = rand(rng,n,1)

  # create primal feasible point -> b
  Strue = generatePosDefMatrix(r,rng)
  strue = vec(Strue)
  b = A*xtrue + strue

  # create dual feasible point -> P,q
  # consider three different types of the quadratic cost
  # if iii <= 15
  #   P = sparse(generatePosDefMatrix(n,rng))
  # elseif iii <=30
  #   P = speye(n)
  # else
    P = spzeros(n,n)
  # end
  Ytrue = generatePosDefMatrix(r,rng)
  ytrue = vec(Ytrue)
  q = -P*xtrue -  A'*ytrue

  # solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, x[1:n])
  @variable(model, S[1:r,1:r],SDP)
  @objective(model, Min, (0.5*x'*P*x+q'*x)[1])
  B = reshape(b,r,r)
  @expression(model, shared, 0)
  for jjj=1:n
    Ai = reshape(A[:,jjj],r,r)
    @expression(model, shared + Ai*x[jjj])
  end
  @constraint(model, shared+S.== B)
  status = JuMP.solve(model)
  println("Objective value: ", getobjectivevalue(model))
  println("x = ", getvalue(x))

  objTrue = getobjectivevalue(model)
  solTrue = getvalue(x)

#   fn = "SDPQuad$(iii).jld"
#   JLD.save(dirPath*fn,"n",n,"m",m,"A",A,"b",b,"P",P,"q",q,"objTrue",objTrue,"solTrue",solTrue)
#   println("$(iii)/$(nn) completed!")
# end


