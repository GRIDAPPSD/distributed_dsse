#!/usr/bin/env julia

using ADNLPModels, NLPModelsIpopt
using CSV


function parse_measurements()
  nodename_nodeidx_map = Dict()
  nnode = 0
  for row in CSV.File("test_4/nodelist.csv", header=false)
    nnode += 1
    nodename_nodeidx_map[row[1]] = nnode
  end
  println("    Total number of nodes: $(nnode)")

  # figure out how many total measurements we have for allocating arrays
  nmeas = (open("test_4/measurements.csv") |> readlines |> length) -1
  println("    Total number of measurements: $(nmeas)")

  vi_measidxs = Vector{Int64}()
  Ti_measidxs = Vector{Int64}()
  measidx_nodeidx_map = Dict()
  rmat = Array{Float64}(undef, nmeas)

  imeas = 0
  for row in CSV.File("test_4/measurements.csv")
    # columns: sensor_type[1],sensor_name[2],node1[3],node2[4],value[5],sigma[6],is_pseudo[7],nom_value[8]

    supported = false
    if row[1] == "vi"
      supported = true
      imeas += 1
      append!(vi_measidxs, imeas)
    elseif row[1] == "Ti"
      supported = true
      imeas += 1
      append!(Ti_measidxs, imeas)
    end

    if supported
      rmat[imeas] = row[6]^2
      measidx_nodeidx_map[imeas] = nodename_nodeidx_map[row[3]]
    end
  end
  println("    Number of vi measurements: $(length(vi_measidxs))")
  println("    Number of Ti measurements: $(length(Ti_measidxs))")

  return nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs
end


# Main

println("Start parsing files...")

nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs = parse_measurements()

println("Done parsing files, start defining optimization problem...")

# just hardwire initial x starting point, x0, for now
x0 = Array{Float64}(undef, 2*nnode)
for i = 1:nnode
  x0[i] = 1.0
end
for i = nnode+1:2*nnode
  x0[i] = 0.0
end

# initialize z here then set values later when iterating over
# measurement_data.csv
zvec = Array{Float64}(undef, nmeas)

vi(x) = sum((zvec[i] - x[measidx_nodeidx_map[i]])^2/rmat[i] for i in vi_measidxs)
Ti(x) = sum((zvec[i] - x[measidx_nodeidx_map[i]+nnode])^2/rmat[i] for i in Ti_measidxs)

f(x) = (length(vi_measidxs)>0 ? vi(x) : 0) + (length(Ti_measidxs)>0 ? Ti(x) : 0)

nlp = ADNLPModel(f, x0)

println("Done with defining optimization problem, start solving it...")

# process each timestamp getting measurement data and calling solver

for row in CSV.File("test_4/measurement_data.csv")
  println("\n================================================================================\n")

  # This logic assumes that the order of measurements (columns in a row)
  # in meaurement_data.csv is the same as measurement order (rows)
  # in measurements.csv
  # If that's not the case, then this will require some fanciness with
  # dictionaries to assign zvec values to the appropriate index
  for imeas in 1:nmeas
    zvec[imeas] = row[imeas+1]
  end
  println("*** Timestamp $(row[1]), using measurements: $(zvec)\n")

  stats = ipopt(nlp)
  print(stats)
end

