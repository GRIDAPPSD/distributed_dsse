#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays


function get_input(zonenum)
  nodename_nodeidx_map = Dict()

  inode = 0
  for row in CSV.File(string("mase_files/nodelist.csv.", zonenum), header=false)
    inode += 1
    nodename_nodeidx_map[row[1]] = inode # this one is needed below
  end
  println("    Total number of nodes: $(inode)")

  measidxs = Dict()
  measidx_nodeidx_map = Dict()
  rmat = Vector{Float64}()

  imeas = 0
  for row in CSV.File(string("mase_files/measurements.csv.", zonenum))
    # columns: sensor_type[1],sensor_name[2],node1[3],node2[4],value[5],sigma[6],is_pseudo[7],nom_value[8]

    stype = row[1]
    if !(stype in keys(measidxs))
      measidxs[stype] = Vector{Int64}()
    end
    imeas += 1
    append!(measidxs[stype], imeas)
    measidx_nodeidx_map[imeas] = nodename_nodeidx_map[row[3]]
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

  println("    measurement node index map: $(measidx_nodeidx_map)")

  # TODO: objective function works on dictionaries rather than potentially
  # faster Julia SparseArray data types. Commented out code below populates
  # YbusG and YbusB SparseArrays.  Tried to plug this into objective function
  # using the nzrange() function in place of keys() for a dictionary, but gave
  # an out of bounds array access error so went back to dictionary
  Ybus = Dict()
  #YbusG = spzeros(Float64, inode, inode)
  #YbusB = spzeros(Float64, inode, inode)
  ibus = 0
  for row in CSV.File(string("mase_files/ysparse.csv.", zonenum))
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
  for row in CSV.File(string("mase_files/vnom.csv.", zonenum))
    if row[1] in keys(nodename_nodeidx_map)
      Vnom[nodename_nodeidx_map[row[1]]] = (row[2], row[3])
      inom += 1
    end
  end
  println("    Vnom number of elements: $(inom)")

  Source = Vector{Int64}()
  for row in CSV.File(string("mase_files/sourcebus.csv.", zonenum), header=false)
    if row[1] in keys(nodename_nodeidx_map)
      append!(Source, nodename_nodeidx_map[row[1]])
    end
  end
  println("    source bus indices: $(Source)")

  measdata = CSV.File(string("mase_files/measurement_data.csv.", zonenum))

  return measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, Source, measdata
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
#  - Source: source bus node indices as a vector
#  - measdata: streaming measurement data for populating zvec

function setup_estimate(measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, Source)
  nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-8,"acceptable_tol"=>1e-8,"max_iter"=>5000,"linear_solver"=>"mumps"))

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
      if i in Source
        @NLconstraint(nlp, v[i] == Vnom[i][1])
      end
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
      if i in Source
        @NLconstraint(nlp, T[i] == deg2rad(start))
      else
        # if this angle constraint is broken up into two NLconstraints, one
        # <= and one >=, bad things will happen with JuMP including crashes
        @NLconstraint(nlp, deg2rad(start-90.0) <= T[i] <= deg2rad(start+90.0))
      end
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


function estimate(nlp, zvec, v, T, measurement)
  # populate zvec with measurement data for given timestep and call solver

  println("\n================================================================================\n")
  # This logic assumes that the order of measurement data (columns in a row)
  # is the same as measurement order (rows) in measurements.csv
  # If that's not the case, then this will require some fanciness with
  # dictionaries to assign zvec values to the appropriate index
  nmeas = length(measurement) - 1
  for imeas in 1:nmeas
    set_value(zvec[imeas], measurement[imeas+1])
  end
  println("*** Timestamp with measurements: $(measurement)\n")

  @time optimize!(nlp)
  solution_summary(nlp, verbose=true)
  println("v = $(value.(v))")
  println("T = $(value.(T))")
end


# Main

println("Start parsing input files...")

measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, Source, measdata = get_input(0)

println("Done parsing input files, start defining optimization problem...")

#nlp, zvec, v, T = setup_estimate(measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom, Source)

#println("Done with defining optimization problem, start solving it...")

#for measurement in measdata
#  estimate(nlp, zvec, v, T, measurement)
#end

