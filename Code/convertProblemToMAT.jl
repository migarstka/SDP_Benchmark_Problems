# Convert all problems from JLD format into MATLAB's MAT format

using FileIO, MAT
fromPath = "../DataFiles/JLD2/QP/QP-Lasso/"
toPath = "../DataFiles/JLD2/QP/QP-Lasso-MATLAB/"

problems = []
  for f in filter(x -> endswith(x, ".jld2"), readdir(fromPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end

  for p in problems
    # load data from JLD
    data = load("$(fromPath)"*"$(p).jld2")

    mfile = matopen(toPath*"/"*p*".mat", "w")
    # iterate over all keys available
    for k in collect(keys(data))
      write(mfile,k,data[k])
    end
    close(mfile)
    println(" "^8*"$(p) completed!")

  end


