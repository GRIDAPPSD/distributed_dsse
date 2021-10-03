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
    # must construct full Ybus, not just lower diagonal elements
    Ybus[row[1]][row[2]] = Ybus[row[2]][row[1]] = complex(row[3], row[4])
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
  x0[i] = 1.0 # initial magnitudes
end
for i = nnode+1:2*nnode
  x0[i] = 0.0 # initial angles
end

# initialize z here then set values later when iterating over
# measurement_data.csv
zvec = Array{Float64}(undef, nz)

# Decomposition primitives for Vi/Ti measurements
vi(x,i) = x[i]
vzi(x,zidx) = vi(x,zidx_nodeidx_map[zidx])

Ti(x,i) = x[i+nnode]
Tzi(x,zidx) = Ti(x,zidx_nodeidx_map[zidx])

# Decomposed objective functions for Vi/Ti measurements
Vi(x) = sum((zvec[zidx] - vzi(x,zidx))^2/rmat[zidx] for zidx in vi_zidxs)
Ti(x) = sum((zvec[zidx] - Tzi(x,zidx))^2/rmat[zidx] for zidx in Ti_zidxs)

# Full form objective functions for Vi/Ti measurements
# Comment either the ones below out or the ones just above
#Vi(x) = sum((zvec[zidx] - x[zidx])^2/rmat[zidx] for zidx in vi_zidxs)
#Ti(x) = sum((zvec[zidx] - x[zidx+nnode])^2/rmat[zidx] for zidx in Ti_zidxs)

# For Pi measurement, from the reference book, "Power System State Estimation 
# Theory and Implementation":
# Pi = vi * sum(vj*(G*cos(T) + B*sin(T)) for j in keys(Ybus[zidx_nodeidx_map[zidx]])
# where
#  vi = x[zidx_nodeidx_map[zidx]]
#  vj = x[j]
#  T = x[zidx_nodeidx_map[zidx] + nnode] - x[j + nnode]
#  G = real(Ybus[zidx_nodeidx_map[zidx]][j])
#  B = imag(Ybus[zidx_nodeidx_map[zidx]][j])

# Decomposition primitives for Pi/Qi measurements
Gij(i,j) = real(Ybus[i][j])
Bij(i,j) = imag(Ybus[i][j])

Gzij(zidx,j) = Gij(zidx_nodeidx_map[zidx],j)
Bzij(zidx,j) = Bij(zidx_nodeidx_map[zidx],j)

vj(x,j) = vi(x,j)

Tij(x,i,j) = Ti(x,i) - Ti(x,j)
Tzij(x,zidx,j) = Tij(x,zidx_nodeidx_map[zidx],j)

# Decomposed objective functions for Pi/Qi measurements
h_Pi(x,zidx) = vzi(x,zidx) * sum(vj(x,j)*(Gzij(zidx,j)*cos(Tzij(x,zidx,j)) + Bzij(zidx,j)*sin(Tzij(x,zidx,j))) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))
Pi(x) = sum((zvec[zidx] - h_Pi(x,zidx))^2/rmat[zidx] for zidx in Pi_zidxs)

h_Qi(x,zidx) = vzi(x,zidx) * sum(vj(x,j)*(Gzij(zidx,j)*sin(Tzij(x,zidx,j)) - Bzij(zidx,j)*cos(Tzij(x,zidx,j))) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))
Qi(x) = sum((zvec[zidx] - h_Qi(x,zidx))^2/rmat[zidx] for zidx in Qi_zidxs)

Pi(x) = sum((zvec[zidx] - h_Pi(x,zidx))^2/rmat[zidx] for zidx in Pi_zidxs)
Qi(x) = sum((zvec[zidx] - h_Qi(x,zidx))^2/rmat[zidx] for zidx in Qi_zidxs)

# Full form objective functions for Pi/Qi measurements
# Comment either the ones below out or the ones just above
#Pi(x) = sum((zvec[zidx] - (x[zidx_nodeidx_map[zidx]] * sum(x[j]*(real(Ybus[zidx_nodeidx_map[zidx]][j])*cos(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode]) + imag(Ybus[zidx_nodeidx_map[zidx]][j])*sin(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))))^2/rmat[zidx] for zidx in Pi_zidxs)
#Qi(x) = sum((zvec[zidx] - (x[zidx_nodeidx_map[zidx]] * sum(x[j]*(real(Ybus[zidx_nodeidx_map[zidx]][j])*sin(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode]) - imag(Ybus[zidx_nodeidx_map[zidx]][j])*cos(x[zidx_nodeidx_map[zidx]+nnode] - x[j+nnode])) for j in keys(Ybus[zidx_nodeidx_map[zidx]]))))^2/rmat[zidx] for zidx in Qi_zidxs)

# Final objective function over all measurement types
f(x) = (length(vi_zidxs)>0 ? Vi(x) : 0) + (length(Ti_zidxs)>0 ? Ti(x) : 0) + (length(Pi_zidxs)>0 ? Pi(x) : 0) + (length(Qi_zidxs)>0 ? Qi(x) : 0)

nlp = ADNLPModel(f, x0)

println("Done with defining optimization problem, start solving it...")

# process each timestamp getting measurement data and calling solver

# MATLAB FPI solver solution for the 4-node test case 
target_solution = [1.0000, 0.971576218875089, 0.956643047805180, 0.950535446461549, 0.0, 0.0, 0.0, 0.0] # FPI solver 10,000 iterations

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

  # Note this needs to be commented out for other than the 4-node test case
  diff_solution = stats.solution - target_solution
  println("\nTarget solution difference:  $(diff_solution)")

  #pisolution = Pi(stats.solution)
  #println("\nPi(solution):  $(pisolution)")

  #for zidx in Pi_zidxs
  #  hi = h_Pi(stats.solution,zidx)
  #  println("\nh_$(zidx)(x*): $(hi)")

  #  hi_exp = h_Pi(target_solution,zidx)
  #  println("h_$(zidx)(x_exp): $(hi_exp)")

  #  for j in keys(Ybus[zidx_nodeidx_map[zidx]])
  #    branchP = vzi(target_solution,zidx) * vj(target_solution,j)*(Gzij(zidx,j)*cos(Tzij(target_solution,zidx,j)) + Bzij(zidx,j)*sin(Tzij(target_solution,zidx,j)))
  #    i = zidx_nodeidx_map[zidx]
  #    println("\tbranch_$(i)_$(j): $(branchP)")
  #  end
  #end
  
end

