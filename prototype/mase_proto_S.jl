#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays


function get_input(zone, shared_nodenames)
  println("    Reading input files for zone: $(zone)")
  nodename_nodeidx_map = Dict()
  nodename = Vector{String}()

  inode = 0
  for row in CSV.File(string("mase_files/nodelist.csv.", zone), header=false)
    inode += 1
    nodename_nodeidx_map[row[1]] = inode
    push!(nodename, row[1])

    if row[1] in keys(shared_nodenames)
      push!(shared_nodenames[row[1]], (zone, inode))
    end
  end
  println("    Total number of nodes: $(inode)")

  measidxs = Dict()
  measidx_nodeidx_map = Dict()
  shared_nodeidx_measidx_map = Dict()
  rmat = Vector{Float64}()

  imeas = 0
  for row in CSV.File(string("mase_files/measurements.csv.", zone))
    # columns: sensor_type[1],sensor_name[2],node1[3],node2[4],value[5],sigma[6],is_pseudo[7],nom_value[8]

    stype = row[1]
    if !(stype in keys(measidxs))
      measidxs[stype] = Vector{Int64}()
    end
    imeas += 1
    append!(measidxs[stype], imeas)

    measidx_nodeidx_map[imeas] = nodename_nodeidx_map[row[3]]

    if row[3] in keys(shared_nodenames) && stype!="Ti"
      nodeidx = nodename_nodeidx_map[row[3]]
      if !(nodeidx in keys(shared_nodeidx_measidx_map))
        shared_nodeidx_measidx_map[nodeidx] = Vector{Int64}()
      end

      append!(shared_nodeidx_measidx_map[nodeidx], imeas)
    end

    append!(rmat, row[6]^2)
  end

  println("    Total number of measurements: $(imeas)")
  if "vi" in keys(measidxs)
    println("    Number of vi measurements: $(length(measidxs["vi"]))")
    println("    vi measurement indices: $(measidxs["vi"])")
  end
  if "Ti" in keys(measidxs)
    println("    Number of Ti measurements: $(length(measidxs["Ti"]))")
    println("    Ti measurement indices: $(measidxs["Ti"])")
  end
  if "Pi" in keys(measidxs)
    println("    Number of Pi measurements: $(length(measidxs["Pi"]))")
    println("    Pi measurement indices: $(measidxs["Pi"])")
  end
  if "Qi" in keys(measidxs)
    println("    Number of Qi measurements: $(length(measidxs["Qi"]))")
    println("    Qi measurement indices: $(measidxs["Qi"])")
  end

  println("    Measurement node index map: $(measidx_nodeidx_map)")
  println("    Shared node measurement index map: $(shared_nodeidx_measidx_map)")

  # TODO: objective function works on dictionaries rather than potentially
  # faster Julia SparseArray data types. Commented out code below populates
  # YbusG and YbusB SparseArrays.  Tried to plug this into objective function
  # using the nzrange() function in place of keys() for a dictionary, but gave
  # an out of bounds array access error so went back to dictionary
  Ybus = Dict()
  #YbusG = spzeros(Float64, inode, inode)
  #YbusB = spzeros(Float64, inode, inode)
  ibus = 0
  for row in CSV.File(string("mase_files/ysparse.csv.", zone))
    if !haskey(Ybus, row[1])
      Ybus[row[1]] = Dict()
    end
    # must construct full Ybus, not just lower diagonal elements
    Ybus[row[1]][row[2]] = Ybus[row[2]][row[1]] = complex(row[3], row[4])
    #YbusG[row[1],row[2]] = YbusG[row[2],row[1]] = row[3]
    #YbusB[row[1],row[2]] = YbusB[row[2],row[1]] = row[4]
    ibus += 1
  end
  println("    Ybus number of lower diagonal elements: $(ibus)")

  Vnom = Dict()
  inom = 0
  for row in CSV.File(string("mase_files/vnom.csv.", zone))
    if row[1] in keys(nodename_nodeidx_map)
      Vnom[nodename_nodeidx_map[row[1]]] = (row[2], row[3])
      inom += 1
    end
  end
  println("    Vnom number of elements: $(inom)")

  source_nodeidxs = Vector{Int64}()
  for row in CSV.File(string("mase_files/sourcebus.csv.", zone), header=false)
    if row[1] in keys(nodename_nodeidx_map)
      append!(source_nodeidxs, nodename_nodeidx_map[row[1]])
    end
  end
  println("    Source node indices: $(source_nodeidxs)\n")

  measdata = CSV.File(string("mase_files/measurement_data.csv.", zone))

  return measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, source_nodeidxs, nodename, shared_nodeidx_measidx_map, measdata
end


# Input:
#  - measidxs: Measurement indices for each type of measurement (vi,Ti,Pi,Qi),
#    dictionary indexed by measurement type with value as vector of measurement
#    indices for the given type
#  - measidx_nodeidx_map: Measurement index number to node index number map
#  - rmat: R measurement covariance matrix, diagonal terms only, square of
#    standard deviation of error of corresponding measurement
#  - Ybus: sparse admittance matrix, 2-dimensional complex value dictionary
#    indexed by node numbers with both lower and upper diagonal values populated
#  - Vnom: nominal voltages, dictionary indexed by node numbers with each
#    element as a tuple with magnitude first and angle (degrees) second
#  - source_nodeidxs: source bus node indices as a vector
#  - measdata: streaming measurement data for populating zvec

# State as of Feb 28:
# * Julia crashes are mostly associated with zones 2 and 4, but not always
# With physical unit inputs:
#   * With the source bus angle equality constraint Julia crashes with
#     tolerance tighther than 1e-8
#   * With only the +/-90 degree inequality angle constraint, it will optimize
#     with any tolerance value
# With per-unit inputs:
#   * objective function values always a couple of orders of magnitude
#     higher than with physical units
#   * With the source bus angle equality constraint Julia crashes with
#     tolerance tighther than 1e-10
#   * With only the +/-90 degree inequality angle constraint, it will optimize
#     with any tolerance value

function setup_estimate(measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, source_nodeidxs)
  nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-8,"acceptable_tol"=>1e-8,"max_iter"=>1000,"linear_solver"=>"mumps"))
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-8))
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer)) # tol default to 1e-8 and acceptable_tol defaults to 1e-6
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-10))
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-12,"acceptable_tol"=>1e-12,"max_iter"=>2000)) # not recommended, lots of iterating with no better results

  nnode = length(Vnom) # get number of nodes from # of Vnom elements
  nmeas = length(rmat) # get number of measurements from # of rmat elements

  # declare z here then set values later during estimate step
  @NLparameter(nlp,zvec[i=1:nmeas] == 0.0)

  # Starting conditions and constraints

  # After investigating what is needed for starting values and constraints,
  # the conclusion is that if we can limit the model problem to feeders where
  # the substation transformer isn't included, we should be able to avoid
  # topological analysis.  But, for the general case we assume we'll have/need
  # all nominal voltages and that is how it is coded below for the prototype.

  # For magnitudes, although we produce as good of a solution without knowing
  # all the nominal voltages and just some of them, having to know any of them
  # means we might as well utilize all of them to give the best shot at
  # finding the correct solution for any model

  @variable(nlp,v[1:nnode])
  for i = 1:nnode
    if i in keys(Vnom)
      set_start_value.(v[i], Vnom[i][1])
# TODO: do we want this source bus magnitude equality constraint?
#      if i in source_nodeidxs
#        @NLconstraint(nlp, v[i] == Vnom[i][1])
#      end
    end
  end

  # Similar to magnitudes, for angles if we need the nominal values at all,
  # we should really use all of them to give us the best shot at finding the
  # correct solution for any model

  @variable(nlp,T[1:nnode])
  for i = 1:nnode
    if i in keys(Vnom)
      start = Vnom[i][2]
      set_start_value.(T[i], deg2rad(start))
# TODO: do we want this source bus angle equality constraint?
#      if i in source_nodeidxs
#        @NLconstraint(nlp, T[i] == deg2rad(start))
#      else
#        # if this angle constraint is broken up into two NLconstraints, one
#        # <= and one >=, bad things will happen with JuMP including crashes
#        @NLconstraint(nlp, deg2rad(start-90.0) <= T[i] <= deg2rad(start+90.0))
#     end
      @NLconstraint(nlp, deg2rad(start-90.0) <= T[i] <= deg2rad(start+90.0))
    end
  end

  # Objective function formulation
  @NLexpression(nlp, vzi[i=1:nmeas], v[measidx_nodeidx_map[i]])

  # vi measurements
  if "vi" in keys(measidxs)
    @NLexpression(nlp, Visum, sum((zvec[i] - vzi[i])^2/rmat[i] for i in measidxs["vi"]))
  else
    @NLexpression(nlp, Visum, 0.0)
  end

  @NLexpression(nlp, Tzi[i=1:nmeas], T[measidx_nodeidx_map[i]])

  # Ti measurements
  if "Ti" in keys(measidxs)
    @NLexpression(nlp, Tisum, sum((zvec[i] - Tzi[i])^2/rmat[i] for i in measidxs["Ti"]))
  else
    @NLexpression(nlp, Tisum, 0.0)
  end

  # common terms for Pi and Qi measurements
  if "Pi" in keys(measidxs) || "Qi" in keys(measidxs)
    @NLexpression(nlp, Tzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], Tzi[i] - T[j])
    @NLexpression(nlp, Gzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], real.(Ybus[measidx_nodeidx_map[i]][j]))
    @NLexpression(nlp, Bzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], imag.(Ybus[measidx_nodeidx_map[i]][j]))
  end

  # Pi measurements
  if "Pi" in keys(measidxs)
    @NLexpression(nlp, h_Pi[i in measidxs["Pi"]], vzi[i] * sum(v[j]*(Gzij[i,j]*cos(Tzij[i,j]) + Bzij[i,j]*sin(Tzij[i,j])) for j in keys(Ybus[measidx_nodeidx_map[i]])))
    @NLexpression(nlp, Pisum, sum((zvec[i] - h_Pi[i])^2/rmat[i] for i in measidxs["Pi"]))
  else
    @NLexpression(nlp, Pisum, 0.0)
  end

  # Qi measurements
  if "Qi" in keys(measidxs)
    @NLexpression(nlp, h_Qi[i in measidxs["Qi"]], vzi[i] * sum(v[j]*(Gzij[i,j]*sin(Tzij[i,j]) - Bzij[i,j]*cos(Tzij[i,j])) for j in keys(Ybus[measidx_nodeidx_map[i]])))
    @NLexpression(nlp, Qisum, sum((zvec[i] - h_Qi[i])^2/rmat[i] for i in measidxs["Qi"]))
  else
    @NLexpression(nlp, Qisum, 0.0)
  end

  # final objective function is sum of components
  @NLobjective(nlp, Min, Visum + Tisum + Pisum + Qisum)

  return nlp, zvec, v, T
end


function estimate(nlp, v, T, nodename, zone)
  # populate zvec with measurement data for given timestep and call solver

  @time optimize!(nlp)
  solution_summary(nlp, verbose=true)
  println("v = $(value.(v))")
  println("T = $(value.(T))")

  println("\n*** Solution for zone $(zone):")
  inode = 0
  toterr = 0.0
  for row in CSV.File(string("mase_files/result_data.csv.", zone), header=false)
    inode += 1
    expected = row[2]
    solution = value.(v[inode])
    pererr = 100.0 * abs(solution - expected)/expected
    toterr += pererr
    println("$(nodename[inode]) v exp: $(expected), sol: $(solution), %err: $(pererr)")
  end
  avgerr = toterr/inode
  println("*** Average %err zone $(zone): $(avgerr)")

  #nnode = length(nodename)
  #for inode in 1:nnode
  #  println("$(nodename[inode]) v exp: XXX, sol: $(value.(v[inode])), %err: XXX")
  #end
end


# Main

println("Start parsing input files...")

shared_nodenames = Dict()
for row in CSV.File("mase_files/sharednodelist.csv", header=false)
  shared_nodenames[row[1]] = []
end

measidxs = Dict()
measidx_nodeidx_map = Dict()
rmat = Dict()
Ybus = Dict()
Vnom = Dict()
source_nodeidxs = Dict()
nodename = Dict()
shared_nodeidx_measidx_map = Dict()
measdata = Dict()

for zone = 0:5
  measidxs[zone], measidx_nodeidx_map[zone], rmat[zone], Ybus[zone], Vnom[zone], source_nodeidxs[zone], nodename[zone], shared_nodeidx_measidx_map[zone], measdata[zone] = get_input(zone, shared_nodenames)
end

Sharednodes = Dict()
Sharedmeas = Dict()
for (key, value) in shared_nodenames
  zone = value[1][1]
  inode = value[1][2]
  if !haskey(Sharednodes, zone)
    Sharednodes[zone] = Dict()
    Sharedmeas[zone] = Dict()
  end
  Sharednodes[zone][inode] = value[2]

  for imeas in shared_nodeidx_measidx_map[zone][inode]
    Sharedmeas[zone][imeas] = value[2]
  end

  zone = value[2][1]
  inode = value[2][2]
  if !haskey(Sharednodes, zone)
    Sharednodes[zone] = Dict()
    Sharedmeas[zone] = Dict()
  end
  Sharednodes[zone][inode] = value[1]

  for imeas in shared_nodeidx_measidx_map[zone][inode]
    Sharedmeas[zone][imeas] = value[1]
  end
end
println("Shared nodenames dictionary: $(shared_nodenames)\n")
println("Sharednodes dictionary: $(Sharednodes)\n")
println("Sharedmeas dictionary: $(Sharedmeas)\n")

println("Done parsing input files, start defining optimization problem...")

nlp = Dict()
zvec = Dict()
v = Dict()
T = Dict()

for zone = 0:5
  nlp[zone], zvec[zone], v[zone], T[zone] = setup_estimate(measidxs[zone], measidx_nodeidx_map[zone], rmat[zone], Ybus[zone], Vnom[zone], source_nodeidxs[zone])
end

println("\nDone defining optimization problem, start solving it...")

# assume all measurement_data files contain the same number of rows/timestamps
nrows = length(measdata[0])
println("number of timestamps to process: $(nrows)")
nnode = length(Vnom[0])

for row = 1:1 # first timestamp only
#for row = 1:nrows # all timestamps

  for zone = 0:5
  #for zone = 2:2
#     if zone!=2 && zone!=4
     # This logic assumes that the order of measurement data (columns in a row)
     # is the same as measurement order (rows) in measurements.csv
     # If that's not the case, then this will require rework
     measurement = measdata[zone][row]
     for imeas in 1:length(measurement)-1
       set_value(zvec[zone][imeas], measurement[imeas+1])
     end
     #println("*** Timestamp with measurements: $(measurement)\n")
#     end
  end

  # first optimization for each zone using vnom data for starting point
  for zone = 0:5
  #for zone = 2:2
#    if zone!=2 && zone!=4
      println("\n================================================================================")
      println("1st optimization for timestamp #$(row), zone: $(zone)\n")
      estimate(nlp[zone], v[zone], T[zone], nodename[zone], zone)
#    end
  end

  # update starting values with v/T solution values from first optimization
  for zone = 0:5
#    if zone!=2 && zone!=4
      set_start_value.(v[zone], value.(v[zone]))
      #set_start_value.(T[zone], value.(T[zone]))
#    end
  end

  # DATA EXCHANGE
  # for Pi and Qi measurement data exchange, need to compute:
  #   S = V.(Ybus x V)*
  #   where S is complex and the real component is the corresponding Pi value
  #   and the imaginary component is the corresponding Qi value
  #   V is the complex solution vector from the first state estimate where the
  #   real component values are the v magnitude values and the imaginary
  #   component values are the T angle values
  #   The "." operation is element-wise multiplication
  #   The "x" operation is regular matrix multiplication
  #   The "*" operation is complex conjugate
  #   Ybus will be created as a dense matrix for the initial implementation


#=
  # exchange shared node values updating v/T starting values
  for (zonedest, Zonenodes) in Sharednodes
    for (nodedest, source) in Zonenodes
      println("starting value sharing destination zone: $(zonedest), node: $(nodedest), source zone: $(source[1]), node: $(source[2]), value: $(value.(v[source[1]][source[2]]))")
      set_start_value.(v[zonedest][nodedest], value.(v[source[1]][source[2]]))
      # don't exchange angle values
      #set_start_value.(T[zonedest][nodedest], value.(T[source[1]][source[2]]))
    end
  end

  # exchange shared node values updating zvec measurement values
  for (zonedest, Zonemeas) in Sharedmeas
    for (measdest, source) in Zonemeas
      println("measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), value: $(value.(v[source[1]][source[2]]))")
      set_value(zvec[zonedest][measdest], value.(v[source[1]][source[2]]))
    end
  end

  # second optimization using shared node values
  for zone = 0:5
#    if zone!=2 && zone!=4
      println("\n================================================================================")
      println("2nd optimization for timestamp #$(row), zone: $(zone)\n")
      estimate(nlp[zone], v[zone], T[zone], nodename[zone], zone)
#    end
  end

  # reset starting values back to Vnom so there is no "timestamp memory"
  # no need to set constraints because those were never updated with
  # shared node data exchange
  for zone = 0:5
#    if zone!=2 && zone!=4
      for i = 1:nnode
        if i in keys(Vnom[zone])
          set_start_value.(v[zone][i], Vnom[zone][i][1])
          set_start_value.(T[zone][i], deg2rad(Vnom[zone][i][2]))
        end
#      end
    end
  end
=#

end

