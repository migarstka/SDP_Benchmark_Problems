
workspace()
include("../../../src/Solver.jl")
include("SparseSDPs.jl")
using OSSDP,SparseSDPs
using Base.Test
rng = MersenneTwister(3232)


function inCone(s,K,tol)

  b = 1
  if K.f > 0
    e = b + K.f - 1
    if maximum(abs.(s[b:e])) > tol
      println("s component not in {0} cone.")
      return false
    end
    b = e + 1
  end

  if K.l > 0
    e = b + K.l - 1
    if minimum(s[K.f+1:K.f+K.l-1]) < -tol
      println("s component not in R+ cone.")
      return false
    end
    b = e+1
  end

  if length(K.q) > 0
    for iii = 1:length(K.q)
      e = b + K.q[iii] - 1
      s_sub = s[b:e]
      if norm(s_sub[2:end],2) -s_sub[1] > -tol
        println("s component not in Lorenz cone nr $(iii).")
        return false
      end
      b = e+1
    end
  end

  b = K.f + K.l + sum(K.q) + 1
  if length(K.s) > 0
    for iii = 1:length(K.s)
      e = b + K.s[iii] - 1
      s_sub = s[b:e]
      cDim = Int(sqrt(K.s[iii]))
      if minimum(eig(reshape(full(s_sub),cDim,cDim))[1]) < -tol
        println("s component not in psd cone nr $(iii).")
        return false
      end
      b = e + 1
    end
  end

  return true
end

# generate a random test problem with banded sdp sparsity pattern
m = 3
K = OSSDPTypes.Cone(2,3,[],[10000 16])
bandwidth = [5 1]
density = 0.45
P,q,A,b,K = generateBandedSDP(rng,m,K,bandwidth,-20.,20.,density);
P = sparse(P)

# one solve to compile functions (to get timings right)
settings = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=3,verbose=false,adaptive_rho=true,scaling=10,decompose=false)
res0,nothing = OSSDP.solve(P,q,A,b,K,settings);


# # solve problem with sparsity turned on
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,eps_abs=1e-5,eps_rel=1e-6,max_iter=2000,verbose=false,adaptive_rho=true,scaling=10,decompose=true)
res_decomposed,ws1,normArr1= OSSDP.solve(P,q,A,b,K,settings1);


# solve problem without sparsity exploitation
settings2 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=2000,verbose=false,adaptive_rho=true,scaling=10,decompose=false)
res,ws2,normArr2 = OSSDP.solve(P,q,A,b,K,settings2);


# # compare results
@testset "Banded Cone comparison" begin
  @test abs(res_decomposed.cost-res.cost)  < 1e-3
  @test inCone(res_decomposed.s,K,1e-5)
  @test inCone(res.s,K,1e-5)
end



