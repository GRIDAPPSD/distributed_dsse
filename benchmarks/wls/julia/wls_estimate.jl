#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays


function get_input()
  nodename_nodeidx_map = Dict()

  inode = 0
  for row in CSV.File("test/nodelist.csv", header=false)
    inode += 1
    nodename_nodeidx_map[row[1]] = inode # this one is needed below
  end
  println("    Total number of nodes: $(inode)")

  measidxs = Dict()
  measidx_nodeidx_map = Dict()
  rmat = Vector{Float64}()

  imeas = 0
  for row in CSV.File("test/measurements.csv")
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

  Ybus = Dict()
  for row in CSV.File("test/ysparse.csv")
    if !haskey(Ybus, row[1])
      Ybus[row[1]] = Dict()
    end
    # must construct full Ybus, not just lower diagonal elements
    Ybus[row[1]][row[2]] = Ybus[row[2]][row[1]] = complex(row[3], row[4])
  end

  Vnom = Dict()
  for row in CSV.File("test/vnom.csv")
    if row[1] in keys(nodename_nodeidx_map)
      Vnom[nodename_nodeidx_map[row[1]]] = (row[2], row[3])
    end
  end

  Source = Vector{Int64}()
  for row in CSV.File("test/sourcebus.csv", header=false)
    if row[1] in keys(nodename_nodeidx_map)
      append!(Source, nodename_nodeidx_map[row[1]])
    end
  end
  println(Source)

  measdata = CSV.File("test/measurement_data.csv")

  return Ybus, Vnom, Source, rmat, measidxs, measidx_nodeidx_map, measdata
end


# Schedule meeting with Fernando for Tuesday 2/8 to show the estimate function
# to start the thinking on how to organize the overall distributed SE
# methodology in terms of agents

# Input:
#  - Ybus
#  - Nominal voltage magnitudes and angles per node, Vnom
#  - Source bus nodes, Source
#  - R covariance matrix, rmat
#  - Measurement indices for each type of measurement (vi, Ti, Pi, Qi), measidxs
#  - Measurement index number to node index number map, measidx_nodeidx_map
#  - Streaming measurement data, measdata, used to populate zvec

function setup_estimate(Ybus, Vnom, Source, rmat, measidxs, measidx_nodeidx_map)
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-14,"acceptable_tol"=>1e-14,"max_iter"=>100000)) # force a whole lot of iterations
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-9,"linear_solver"=>"mumps")) # fast optimization with decent accuracy
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-10,"max_iter"=>10000,"linear_solver"=>"ma27"))
  #nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-10,"max_iter"=>10000,"linear_solver"=>"mumps"))
  nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-8,"max_iter"=>5000,"linear_solver"=>"mumps"))

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


function estimate(nlp, zvec, v, T, measdata)
  # process measurement data for each timestep, calling solver

  for measurement in measdata
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
end


# Main

println("Start parsing input files...")

Ybus, Vnom, Source, rmat, measidxs, measidx_nodeidx_map, measdata = get_input()

println("Done parsing input files, start defining optimization problem...")

nlp, zvec, v, T = setup_estimate(Ybus, Vnom, Source, rmat, measidxs, measidx_nodeidx_map)

println("Done with defining optimization problem, start solving it...")

estimate(nlp, zvec, v, T, measdata)

