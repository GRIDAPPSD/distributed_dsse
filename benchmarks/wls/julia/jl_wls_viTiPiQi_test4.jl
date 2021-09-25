#!/usr/bin/env julia

println("Started, before using statements")
using ADNLPModels, NLPModelsIpopt
using CSV
println("After using statements")


function parse_measurements()
  println("    In parse_measurements function")
  nodename_nodeidx_map = Dict()
  nnode = 0
  for row in CSV.File("test/nodelist.csv", header=false)
    nnode += 1
    nodename_nodeidx_map[row[1]] = nnode
  end
  println("    Total number of nodes: $(nnode)")

  # figure out how many total measurements we have for allocating arrays
  nmeas = (open("test/measurements.csv") |> readlines |> length) -1
  println("    Total number of measurements: $(nmeas)")

  vi_zidxs = Vector{Int64}()
  Ti_zidxs = Vector{Int64}()
  Pi_zidxs = Vector{Int64}()
  Qi_zidxs = Vector{Int64}()
  zidx_nodeidx_map = Dict()
  rmat = Array{Float64}(undef, nmeas)

  zidx = 0
  for row in CSV.File("test/measurements.csv")
    # columns: sensor_type[1],sensor_name[2],node1[3],node2[4],value[5],sigma[6],is_pseudo[7],nom_value[8]

    supported = false
    if row[1] == "vi"
      supported = true
      zidx += 1
      append!(vi_zidxs, zidx)
    elseif row[1] == "Ti"
      supported = true
      zidx += 1
      append!(Ti_zidxs, zidx)
    elseif row[1] == "Pi"
      supported = true
      zidx += 1
      append!(Pi_zidxs, zidx)
    elseif row[1] == "Qi"
      supported = true
      zidx += 1
      append!(Qi_zidxs, zidx)
    end

    if supported
      rmat[zidx] = row[6]^2
      zidx_nodeidx_map[zidx] = nodename_nodeidx_map[row[3]]
    end
  end
  println("    Number of vi measurements: $(length(vi_zidxs))")
  println("    Number of Ti measurements: $(length(Ti_zidxs))")
  println("    Number of Pi measurements: $(length(Pi_zidxs))")
  println("    Number of Qi measurements: $(length(Qi_zidxs))")

  Ybus = Dict()
  for row in CSV.File("test/ysparse.csv")
    if !haskey(Ybus, row[1])
      Ybus[row[1]] = Dict()
    end
    Ybus[row[1]][row[2]] = complex(row[3], row[4])
  end
  #println(Ybus)
  #ccall(:jl_exit, Cvoid, (Int32,), 86)

  return nnode, nmeas, rmat, zidx_nodeidx_map, vi_zidxs, Ti_zidxs, Pi_zidxs, Qi_zidxs, Ybus
end


# Main

println("Start parsing files...")

nnode, nz, rmat, zidx_nodeidx_map, vi_zidxs, Ti_zidxs, Pi_zidxs, Qi_zidxs, Ybus = parse_measurements()

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
zvec = Array{Float64}(undef, nz)

vi(x) = sum((zvec[zidx] - x[zidx_nodeidx_map[zidx]])^2/rmat[zidx] for zidx in vi_zidxs)
Ti(x) = sum((zvec[zidx] - x[zidx_nodeidx_map[zidx]+nnode])^2/rmat[zidx] for zidx in Ti_zidxs)

# C++ SE code based Pi(x) formulation
#i = zidx_nodeidx_map[zidx] 
# Ybus is row and column indexed (dict of dict)
# If using this formulation, need to get rid of Ybus for i=j (diagonal elements)
#sum((vi*vi/(ai*ai) * g) - (vi*vj/(ai*aj) * (g*cos(T) + b*sin(T))) for j in Ybus[i].keys) + vi*vi * gii # iterate over all dict keys

# Power System State Estimation Theory and Implementation book Pi(x) formulation
#vi * sum(vj*(G*cos(T) + B*sin(T)) for j in Ybus[i].keys)
#where
#  vi = x[zidx_nodeidx_map[zidx]]
#  vj = x[j]
#  T = x[zidx_nodeidx_map[zidx] + nnode] - x[j + nnode]
#  G = real(Ybus[zidx_nodeidx_map[zidx]][j])
#  B = imag(Ybus[zidx_nodeidx_map[zidx]][j])
#therefore
  #PiMeasFunc = x[zidx_nodeidx_map[zidx]] * sum(x[j]*(real(Ybus[zidx_nodeidx_map[zidx]][j])*cos(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) + (imag(Ybus[zidx_nodeidx_map[zidx]][j])*sin(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))
  #Pi(x) = sum((zvec[zidx] - PiMeasFunc)^2/rmat[zidx] for zidx in Pi_zidxs)

Pi(x) = sum((zvec[zidx] - (x[zidx_nodeidx_map[zidx]] * sum(x[j]*(real(Ybus[zidx_nodeidx_map[zidx]][j])*cos(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) + (imag(Ybus[zidx_nodeidx_map[zidx]][j])*sin(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))))^2/rmat[zidx] for zidx in Pi_zidxs)

Qi(x) = sum((zvec[zidx] - (x[zidx_nodeidx_map[zidx]] * sum(x[j]*(real(Ybus[zidx_nodeidx_map[zidx]][j])*sin(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) - (imag(Ybus[zidx_nodeidx_map[zidx]][j])*cos(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))))^2/rmat[zidx] for zidx in Pi_zidxs)

f(x) = (length(vi_zidxs)>0 ? vi(x) : 0) + (length(Ti_zidxs)>0 ? Ti(x) : 0) + (length(Pi_zidxs)>0 ? Pi(x) : 0) + (length(Qi_zidxs)>0 ? Qi(x) : 0)

nlp = ADNLPModel(f, x0)

println("Done with defining optimization problem, start solving it...")

# process each timestamp getting measurement data and calling solver

for row in CSV.File("test/measurement_data.csv")
  println("\n================================================================================\n")

  # This logic assumes that the order of measurements (columns in a row)
  # in meaurement_data.csv is the same as measurement order (rows)
  # in measurements.csv
  # If that's not the case, then this will require some fanciness with
  # dictionaries to assign zvec values to the appropriate index
  for zidx in 1:nz
    zvec[zidx] = row[zidx+1]
  end
  println("*** Timestamp $(row[1]), using measurements: $(zvec)\n")

  #stats = ipopt(nlp, print_level=10)
  stats = ipopt(nlp)
  print(stats)
  println("\nFull solution:  $(stats.solution)")
end

