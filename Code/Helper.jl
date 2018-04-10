module Helper

export generatePosDefMatrix

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
end #MODULE