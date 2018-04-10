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

function vecPos(n,i,j)
  ((i > n) || (j > n)) && error("Your index is outside the matrix dimension.")
  return (j-1)*n+i
end

 for iii =1:1:nn
  # number of polynomial variables
  np = 2
  # degree of polynomial
  d = 4
  # randomly create 2d+1 polynomial coefficients starting with p_0 up to p_2d
  p = rand(rng,-5:5,2*d+1,1)
  #p = [5.;4;6;4;1]

  # determine sdp problem matrices
  n = d+1
  P = spzeros(n^2,n^2)
  q = spzeros(n^2,1)
  Aa = spzeros(2*d+1,n^2)
  ba = p

  # A matrix for univariate case
  for k=0:2*d
    for i=0:d, j=0:d
      if i+j == k
        if i<j
          Aa[k+1,vecPos(n,j+1,i+1)] = Aa[k+1,vecPos(n,j+1,i+1)] + 1
        else
          Aa[k+1,vecPos(n,i+1,j+1)] = Aa[k+1,vecPos(n,i+1,j+1)] + 1
        end
      end
    end
  end


  # # A matrix for bivariate case (x,y)
  # for k=0:2*d
  #   for y=0:d, x=0:d
  #     if i+j == k
  #       if i<j
  #         Aa[k+1,vecPos(n,j+1,i+1)] = Aa[k+1,vecPos(n,j+1,i+1)] + 1
  #       else
  #         Aa[k+1,vecPos(n,i+1,j+1)] = Aa[k+1,vecPos(n,i+1,j+1)] + 1
  #       end
  #     end
  #   end
  # end

  #solve accurately once with mosek
  model = Model(solver=MosekSolver())
  @variable(model, Q[1:n,1:n], SDP)
  q = vec(Q)
  @objective(model, Min, 0)
  @constraint(model,Aa*q.==ba)
  # @constraint(model,Q[3,3] == 1)
  # @constraint(model,2*Q[2,3] == 4)
  # @constraint(model,Q[2,2] + 2*Q[1,3] == 6)
  # @constraint(model,2*Q[1,2] == 4)
  # @constraint(model,Q[1,1] == 5)
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


