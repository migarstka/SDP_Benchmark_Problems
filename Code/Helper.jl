module HelperFunctions

export generatePosDefMatrix, transposeVectorized

  # generate a random pos def matrix with eigenvalues between 0.1 and 2
  function generatePosDefMatrix(n::Int64,rng)
      X = rand(rng,n,n)
      # any real square matrix can be QP decomposed into a orthogonal matrix and an uppertriangular matrix R
      Q, R = qr(X)
      eigs = rand(rng,n).*(2.-0.1) .+ 0.1
      X = Q*diagm(eigs)*Q'
      X = 0.5*(X+X')
      return X
  end

  # gives you the position of Aij in vec(A) where A in R^mxn
  function vecPos(m,n,i,j)
    ((i > m) || (j > n)) && error("Your indizes are outside the matrix dimension.")
    return Int((j-1)*m+i)
  end

  # permutes the cols of matrix A, s.t. A*vec(B) == A[:,p]*vec(B'), with B in R^{m_b x n_b}
  function transposeVectorized(A,mb,nb)
    # size of A
    m,n = size(A)

    # find permutation vector p
    k = 0
    p = zeros(Int64,n)
    for j=1:nb, i=1:mb
      k += 1
      p[vecPos(nb,mb,j,i)] = k
    end
    #println("Permutation vector: p=$(p)")
    # permute cols of A
    return A[:,p]

  end

end #MODULE