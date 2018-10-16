# File used to copy the MOSEK and SCS time from the Julia v0.6 compatible JLD2 problem files to the Julia v1.0 compatible files.


using FileIO, JLD2

dest_dir = "../DataFiles/JLD2/SDP/ClosestCorr-Benchmark/"
target_dir = "../DataFiles/JLD2/SDP/ClosestCorr-Benchmark-v1/"

# file assumes same file names
  files = []
  for f in filter(x -> endswith(x, ".jld2"), readdir(dirPath))
      push!(files,String(f))
  end
  nn = length(files)

  for (iii,file) = enumerate(files)
    dest_data = load(dest_dir * file)
    MOSEKcost = dest_data["MOSEKcost"]
    MOSEKtime = dest_data["MOSEKtime"]
    SCScost = dest_data["SCScost"]
    SCStime = dest_data["SCStime"]

    target_data = load(target_dir * file)
    target_data["MOSEKcost"] = MOSEKcost
    target_data["MOSEKtime"] = MOSEKtime
    target_data["SCScost"] = SCScost
    target_data["SCStime"] = SCStime
    target_keys = keys(target_data)

    JLD2.jldopen(target_dir * file, "w") do target_file
      for key in target_keys
        target_file[key] = target_data[key]
      end
    end
    println("$(iii)/$(nn): File: $(file) completed!")
  end