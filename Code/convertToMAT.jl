# Convert all problems from JLD format into MATLAB's MAT format

using JLD, MAT
fromPath = "../DataFiles/Julia/"
toPath = "../DataFiles/MATLAB/"
existingFolders = readdir(toPath)
problemTypes = []
for f in filter(x -> !startswith(x, "."), readdir(fromPath))
    f = split(f,".")[1]
    push!(problemTypes,String(f))
end

for pType in problemTypes
    println(">>> Current ProblemType: $(pType).")
    subDirPath = fromPath * pType *"/"

    problems = []
    for f in filter(x -> endswith(x, ".jld"), readdir(subDirPath))
        f = split(f,".")[1]
        push!(problems,String(f))
    end

    for p in problems
      # load data from JLD
      data = JLD.load("$(subDirPath)"*"$(p).jld")
      P = data["P"]
      q = data["q"]
      A = data["A"]
      b = data["b"]
      m = data["m"]
      n = data["n"]
      objTrue = data["objTrue"]

      # check if folder already exists in MATLAB dir, otherwise create new
      if !in(pType,existingFolders)
        mkdir(toPath*pType*"/")
        push!(existingFolders,pType)
        end
        mfile = matopen(toPath*pType*"/"*p*".mat", "w")
      # iterate over all keys available
      for k in collect(keys(data))
        write(mfile,k,data[k])
      end
      close(mfile)
      println(" "^8*"$(p) completed!")

    end

end
