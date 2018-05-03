# File that loads test problems and intentionally scales the data in a bad way
# Creates 3^4 different scaled versions per test problem
# Scale P, q, A, b with scalars, matrizes from log uniform distribution between
# -1e-3 and 1e3


# workspace()
# include("../../src/Solver.jl")

# load data files
# using JLD

rng = MersenneTwister(12345)

function createRandomSymmetricDiag(rng,n)
  r = Int(sqrt(n))
  d = 10.^(rand(rng,r)*(1.5-(-1.5))-1.5)
  return diag(kron(diagm(d),diagm(d)))
end

# find available problem types in SDP_Benchmark_Problems folder
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
resultPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Scaled/"
existingFolders = readdir(probPath)
problemTypes = []
for f in filter(x -> !startswith(x, "."), readdir(probPath))
    f = split(f,".")[1]
    push!(problemTypes,String(f))
end
# filter!(x->!in(x,["SmallestCircle";"Lovasz";"MIQO"]),problemTypes)
filter!(x->in(x,["SDPandEqConstraints"]),problemTypes)
println(">>> $(length(problemTypes)) Problem Type(s) detected!")

scalingFactor = [1e-3;1;1e3]





# load each problem type
for pType in problemTypes
  println(">>> ProblemType: $(pType)")
  dirPath = probPath*pType*"/"

  # find all file names in problem type folder
  problems = []
  for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end
  # nn = length(problems)
    nn = 1

  # loop over all problems in problem type folder
  for iii =1:1:nn
    pPath = resultPath*pType*"/"
    !ispath(pPath) && mkdir(pPath)



    gc()
    problem = problems[iii]
    data = JLD.load("$(dirPath)"*"$(problem).jld")
    P = data["P"]
    q = data["q"]
    A = data["A"]
    b = data["b"]
    Kf = data["Kf"]
    Kl = data["Kl"]
    Kq = data["Kq"]
    Ks = data["Ks"]

    m,n = size(A)
    kkk = 1
    # scale P,q,A,b
    for ppp=1:3
      Ps = sparse(P*scalingFactor[ppp])
      for qqq=1:3
        qs = q*scalingFactor[qqq]
        for aaa=1:2
          if aaa == 1
            scaling_diag = 10.^(rand(rng,Kf+Kl)*(3-(-3))-3)
            for mmm in Kq
              scaling_diag = vcat(scaling_diag,createRandomSymmetricDiag(rng,mmm))
            end
            for sss in Ks
              scaling_diag = vcat(scaling_diag,createRandomSymmetricDiag(rng,sss))
            end

            As = sparse(diagm(scaling_diag)*A)
          else
            As = A
          end
          for bbb=1:3
            bs = b*scalingFactor[bbb]

            # save scaled data in dictionary
            dataScaled = Dict()
            for k in collect(keys(data))
              if k != "P" && k != "q" && k!= "b" && k !="A" && k != "objTrue" && k != "solTrue"
                dataScaled[k] = data[k]
              end
            end
            dataScaled["P"] = Ps
            dataScaled["q"] = qs
            dataScaled["A"] = As
            dataScaled["b"] = bs


            # calculate optimal solution with MOSEK




            #save file
            nr = "$(kkk)"
            if kkk < 10
              nr = "0$(kkk)"
            end
            fn = problem*"-"*nr*".jld"

            JLD.save(pPath*fn,dataScaled)
            println("$(kkk)/$(nn*3^3*2) of problem $(problem) completed!")
            kkk+=1
          end
        end
      end
    end
  end
end