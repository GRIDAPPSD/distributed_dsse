#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays

test_dir = "mase_files_pu"

function get_input(zone, shared_nodenames)
  println("    Reading input files for zone: $(zone)")
  nodename_nodeidx_map = Dict()
  nodenames = Vector{String}()

  inode = 0
  for row in CSV.File(string(test_dir, "/nodelist.csv.", zone), header=false)
    inode += 1
    nodename_nodeidx_map[row[1]] = inode
    push!(nodenames, row[1])

    if row[1] in keys(shared_nodenames)
      push!(shared_nodenames[row[1]], (zone, inode))
    end
  end
  println("    Total number of nodes: $(inode)")

  # must initialize these 1st and 2nd estimate data structures separately so
  # they don't point to the same thing
  measidxs1 = Dict()
  measidxs2 = Dict()
  measidx1_nodeidx_map = Dict()
  measidx2_nodeidx_map = Dict()
  shared_nodeidx_measidx2_map = Dict()
  rmat1 = Vector{Float64}()
  rmat2 = Vector{Float64}()

  imeas1 = 0
  for row in CSV.File(string(test_dir, "/measurements.csv.", zone))
    # columns: sensor_type[1],sensor_name[2],node1[3],node2[4],value[5],sigma[6],is_pseudo[7],nom_value[8]

    stype = row[1]
    if !(stype in keys(measidxs1))
      measidxs1[stype] = Vector{Int64}()
      measidxs2[stype] = Vector{Int64}()
    end
    imeas1 += 1
    append!(measidxs1[stype], imeas1)
    append!(measidxs2[stype], imeas1)

    measidx1_nodeidx_map[imeas1] = measidx2_nodeidx_map[imeas1] = nodename_nodeidx_map[row[3]]

    append!(rmat1, row[6]^2)
    append!(rmat2, row[6]^2)

    if row[3] in keys(shared_nodenames)
      nodeidx = nodename_nodeidx_map[row[3]]
      if !(nodeidx in keys(shared_nodeidx_measidx2_map))
        shared_nodeidx_measidx2_map[nodeidx] = Vector{Int64}()
      end

      append!(shared_nodeidx_measidx2_map[nodeidx], imeas1)
    end
  end

  imeas2 = imeas1

  # for the 2nd estimates, make sure there are vi, Pi, and Qi measurements
  # for each shared node.  If there aren't add any that are missing
  for (nodeidx, measidx_vec) in shared_nodeidx_measidx2_map
    viFlag = false
    PiFlag = false
    QiFlag = false
    for measidx in measidx_vec
      if measidx in measidxs2["vi"]
        viFlag = true
      elseif measidx in measidxs2["Pi"]
        PiFlag = true
      elseif measidx in measidxs2["Qi"]
        QiFlag = true
      end
    end

    if !viFlag
      println("***missing vi measurement for node index: $(nodeidx)")
      # create new measurement
      imeas2 += 1
      append!(measidxs2["vi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 0.01^2)
    end
    if !PiFlag
      println("***missing Pi measurement for node index: $(nodeidx)")
      # create new measurement
      imeas2 += 1
      append!(measidxs2["Pi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 1.0^2)
    end
    if !QiFlag
      println("***missing Qi measurement for node index: $(nodeidx)")
      # create new measurement
      imeas2 += 1
      append!(measidxs2["Qi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 0.5^2)
    end
  end

  println("    Total number of 1st estimate measurements: $(imeas1)")
  if "vi" in keys(measidxs1)
    println("    Number of vi measurements: $(length(measidxs1["vi"]))")
    println("    vi measurement indices: $(measidxs1["vi"])")
  end
  if "Ti" in keys(measidxs1)
    println("    Number of Ti measurements: $(length(measidxs1["Ti"]))")
    println("    Ti measurement indices: $(measidxs1["Ti"])")
  end
  if "Pi" in keys(measidxs1)
    println("    Number of Pi measurements: $(length(measidxs1["Pi"]))")
    println("    Pi measurement indices: $(measidxs1["Pi"])")
  end
  if "Qi" in keys(measidxs1)
    println("    Number of Qi measurements: $(length(measidxs1["Qi"]))")
    println("    Qi measurement indices: $(measidxs1["Qi"])")
  end

  println("    Total number of 2nd estimate measurements: $(imeas2)")
  if "vi" in keys(measidxs2)
    println("    Number of vi measurements: $(length(measidxs2["vi"]))")
    println("    vi measurement indices: $(measidxs2["vi"])")
  end
  if "Ti" in keys(measidxs2)
    println("    Number of Ti measurements: $(length(measidxs2["Ti"]))")
    println("    Ti measurement indices: $(measidxs2["Ti"])")
  end
  if "Pi" in keys(measidxs2)
    println("    Number of Pi measurements: $(length(measidxs2["Pi"]))")
    println("    Pi measurement indices: $(measidxs2["Pi"])")
  end
  if "Qi" in keys(measidxs2)
    println("    Number of Qi measurements: $(length(measidxs2["Qi"]))")
    println("    Qi measurement indices: $(measidxs2["Qi"])")
  end

  println("    Measurement node index map: $(measidx2_nodeidx_map)")
  println("    Shared node measurement index map: $(shared_nodeidx_measidx2_map)")

  # TODO: objective function works on dictionaries rather than potentially
  # faster Julia SparseArray data types. Commented out code below populates
  # YbusG and YbusB SparseArrays.  Tried to plug this into objective function
  # using the nzrange() function in place of keys() for a dictionary, but gave
  # an out of bounds array access error so went back to dictionary
  Ybus = Dict()
  Ybusp = spzeros(ComplexF64, inode, inode)
  #YbusG = spzeros(Float64, inode, inode)
  #YbusB = spzeros(Float64, inode, inode)
  ibus = 0
  for row in CSV.File(string(test_dir, "/ysparse.csv.", zone))
    if !haskey(Ybus, row[1])
      Ybus[row[1]] = Dict()
    end
    # must construct full Ybus, not just lower diagonal elements
    Ybus[row[1]][row[2]] = Ybus[row[2]][row[1]] = Ybusp[row[1],row[2]] = Ybusp[row[2],row[1]] = complex(row[3], row[4])
    #YbusG[row[1],row[2]] = YbusG[row[2],row[1]] = row[3]
    #YbusB[row[1],row[2]] = YbusB[row[2],row[1]] = row[4]
    ibus += 1
  end
  println("    Ybus number of lower diagonal elements: $(ibus)")

  Vnom = Dict()
  inom = 0
  for row in CSV.File(string(test_dir, "/vnom.csv.", zone))
    if row[1] in keys(nodename_nodeidx_map)
      Vnom[nodename_nodeidx_map[row[1]]] = (row[2], row[3])
      inom += 1
    end
  end
  println("    Vnom number of elements: $(inom)")

  measdata = CSV.File(string(test_dir, "/measurement_data.csv.", zone))

  return measidxs1, measidxs2, measidx1_nodeidx_map, measidx2_nodeidx_map, rmat1, rmat2, Ybus, Ybusp, Vnom, nodenames, nodename_nodeidx_map, shared_nodeidx_measidx2_map, measdata
end


function setup_data_sharing(shared_nodenames, shared_nodeidx_measidx2_map)
  #Sharednodes = Dict()
  Sharedmeas = Dict()
  SharedmeasAlt = Dict()
  for (key, value) in shared_nodenames
    zone = value[1][1]
    inode = value[1][2]
    if !haskey(Sharedmeas, zone)
      #Sharednodes[zone] = Dict()
      Sharedmeas[zone] = Dict()
      SharedmeasAlt[zone] = Dict()
    end
    #Sharednodes[zone][inode] = value[2]

    for imeas in shared_nodeidx_measidx2_map[zone][inode]
      Sharedmeas[zone][imeas] = value[2]
      if zone == 0
        SharedmeasAlt[zone][imeas] = value[2]
      else
        SharedmeasAlt[zone][imeas] = value[1]
      end
    end

    zone = value[2][1]
    inode = value[2][2]
    if !haskey(Sharedmeas, zone)
      #Sharednodes[zone] = Dict()
      Sharedmeas[zone] = Dict()
      SharedmeasAlt[zone] = Dict()
    end
    #Sharednodes[zone][inode] = value[1]

    for imeas in shared_nodeidx_measidx2_map[zone][inode]
      Sharedmeas[zone][imeas] = value[1]
      if zone == 0
        SharedmeasAlt[zone][imeas] = value[1]
      else
        SharedmeasAlt[zone][imeas] = value[2]
      end
    end
  end
  println("Shared nodenames dictionary: $(shared_nodenames)\n")
  #println("Sharednodes dictionary: $(Sharednodes)\n")
  println("Sharedmeas dictionary: $(Sharedmeas)\n")
  println("SharedmeasAlt dictionary: $(SharedmeasAlt)\n")

  return Sharedmeas, SharedmeasAlt
end


function perform_data_sharing(Ybusp, Sharedmeas, SharedmeasAlt, measidxs2, v1, T1, zvec2)
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
  S = Dict()
  for zone = 0:5
    nnode = length(v1[zone]) # get number of nodes from # of v1 elements

    V = Array{ComplexF64}(undef, nnode)
    for inode = 1:nnode
      V[inode] = value.(v1[zone][inode]) * exp(value.(T1[zone][inode])*1im)
    end

    S[zone] = V .* conj(Ybusp[zone] * V)
    println("\nZone: $(zone), nnode: $(nnode), S: $(S[zone])")
  end

  # exchange shared node values updating zvec measurement values
  #for (zonedest, Zonemeas) in SharedmeasAlt
  for (zonedest, Zonemeas) in Sharedmeas
    for (measdest, source) in Zonemeas
      #println("OLD measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), v1 value: $(value.(v1[source[1]][source[2]]))")
      if "vi" in keys(measidxs2[zonedest]) && measdest in measidxs2[zonedest]["vi"]
        println("vi measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), v1 value: $(value.(v1[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], value.(v1[source[1]][source[2]]))
        println("*** zone $(zonedest) zvec2 vi exchange for measurement $(measdest): $(value.(v1[source[1]][source[2]]))")
      elseif "Ti" in keys(measidxs2[zonedest]) && measdest in measidxs2[zonedest]["Ti"]
        println("Ti measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), T1 value: $(value.(T1[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], value.(T1[source[1]][source[2]]))
        println("*** zone $(zonedest) zvec2 Ti exchange for measurement $(measdest): $(value.(T1[source[1]][source[2]]))")
      elseif "Pi" in keys(measidxs2[zonedest]) && measdest in measidxs2[zonedest]["Pi"]
        println("Pi measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), -S real value: $(-1*real(S[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], -1*real(S[source[1]][source[2]]))
        println("*** zone $(zonedest) zvec2 Pi exchange for measurement $(measdest): $(-1*real(S[source[1]][source[2]]))")
      elseif "Qi" in keys(measidxs2[zonedest]) && measdest in measidxs2[zonedest]["Qi"]
        println("Qi measurement value sharing destination zone: $(zonedest), meas: $(measdest), source zone: $(source[1]), node: $(source[2]), -S imag value: $(-1*imag(S[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], -1*imag(S[source[1]][source[2]]))
        println("*** zone $(zonedest) zvec2 Qi exchange for measurement $(measdest): $(-1*imag(S[source[1]][source[2]]))")
      end
    end
  end
end


function buildZonegraph(parzone, shared_nodenames, Zonenodes, Zonegraph)
  println("buildZonegraph parent zone: $(parzone)")
  #fake = 0
  for node in Zonenodes[parzone]
    for (zone, nodeidx) in shared_nodenames[node]
      if zone != parzone && !(zone in keys(Zonegraph))
        if !(parzone in keys(Zonegraph))
          Zonegraph[parzone] = Array{Tuple{Int64,Int64},1}()
        end

        buildZonegraph(zone, shared_nodenames, Zonenodes, Zonegraph)

        println("buildZonegraph Zonegraph[$(parzone)] add child zone: $(zone)")
        if zone in keys(Zonegraph)
          push!(Zonegraph[parzone], (zone, length(Zonegraph[zone])))
        else
          push!(Zonegraph[parzone], (zone, 0))
          #fake += 1
          #TODO given fake counts, see if we can sort on that tuple element
          # or if we need to go to a different data structure besides an array
          # of tuples to allow sorting
          #push!(Zonegraph[parzone], (zone, fake))
        end
      end
    end
  end
end


function buildZoneorder(parzone, Zonegraph, Zoneorder)
  println("buildZoneorder parent zone: $(parzone)")
  if parzone in keys(Zonegraph)
    # order children zones by the count of children each of those have
    order = sort(Zonegraph[parzone], by=x->x[2], rev=true)
    println("buildZoneorder sorted: $(order)")
    # first, add all the children zones in their order
    for (zone, cnt) in order
      append!(Zoneorder, zone)
    end
    # second, traverse those children zones and do the same recursive ordering
    for (zone, cnt) in order
      buildZoneorder(zone, Zonegraph, Zoneorder)
    end
  end
end


function setup_angle_passing(nodename_nodeidx_map, shared_nodenames)
  # first task is determining the zone ordering

  # determine the system reference zone from which zone source nodes are in
  sysref_flag = false
  sysref_zone = 0
  sysref_node = ""
  for row in CSV.File(test_dir*"/sourcenodes.csv", header=false)
    node = row[1]
    for zone = 0:5
      if node in keys(nodename_nodeidx_map[zone])
        if sysref_flag
          println("WARNING: found source nodes in multiple zones")
        else
          sysref_flag = true
          sysref_zone = zone
          sysref_node = node
        end
      end
    end
  end

  if sysref_flag
    println("Found system reference node: $(sysref_node), in zone: $(sysref_zone)")
  else
    println("WARNING: system reference node and zone not found based on source nodes")
  end

  # build a graph of the zones linked to other zones by shared nodes
  # shared_nodenames has the info needed to build this
  println("shared_nodenames: $(shared_nodenames)")

  # for each zone create a list of shared nodes
  Zonenodes = Dict()
  for (node, zonepairs) in shared_nodenames
    for zonepair in zonepairs
      zone = zonepair[1]
      if !(zone in keys(Zonenodes))
        Zonenodes[zone] = Vector{String}()
      end
      push!(Zonenodes[zone], node)
    end
  end
  println("Shared nodes per zone, Zonenodes: $(Zonenodes)")

  Zonegraph = Dict()
  # invoke recursive function to build the graph of connected zones
  # pass the system reference zone and recursion will build the rest
  buildZonegraph(sysref_zone, shared_nodenames, Zonenodes, Zonegraph)
  println("Connected zones graph, Zonegraph: $(Zonegraph)")

  Zoneorder = Vector{Int64}()
  append!(Zoneorder, sysref_zone) # system reference zone is always first
  # traverse the zone graph recursively starting from the system reference
  # zone to build the full zone ordering
  buildZoneorder(sysref_zone, Zonegraph, Zoneorder)
  println("Zone ordering vector, Zoneorder: $(Zoneorder)")

  # create a dictionary to quickly lookup order by zone
  iorder = 0
  ZoneorderDict = Dict()
  # iterate over Zoneorder backwards so the higher priority zones
  # get larger values
  for zone in Iterators.Reverse(Zoneorder)
    iorder += 1
    ZoneorderDict[zone] = iorder
  end
  println("Zone ordering dictionary, ZoneorderDict: $(ZoneorderDict)")

  # second task is determining the zone reference nodes

  Zonerefinfo = Dict()
  for zone = 0:5
    if zone == sysref_zone
      # for the system reference zone, the zone reference node is always
      # the system reference node
      Zonerefinfo[zone] = (sysref_node, nothing, nothing)
    else
      # determine what the shared nodes are for this zone to find the one
      # that is the zone reference node
      # check each shared node to see which has the highest zone order for
      # the other zones where it is shared
      max_priority = 0
      max_node = max_zone = max_idx = nothing
      for node in Zonenodes[zone]
        # find other zone that node is shared with
        for (shared_zone, shared_idx) in shared_nodenames[node]
          if shared_zone!=zone && ZoneorderDict[shared_zone]>max_priority
            max_priority = ZoneorderDict[shared_zone]
            max_node = node
            max_zone = shared_zone
            max_idx = shared_idx
          end
        end
      end
      Zonerefinfo[zone] = (max_node, max_zone, max_idx)
    end
  end
  println("Zone reference info, Zonerefinfo: $(Zonerefinfo)")

  return Zoneorder, Zonerefinfo
end


function perform_angle_passing(T2, Zoneorder, Zonerefinfo, nodename_nodeidx_map, nodenames)
  # store the updated angles after reference angle passing in a new data
  # structure because I can't update the JuMP T2 solution vector
  T2_updated = Dict()

  last_ref_angle = 0.0
  for zone in Zoneorder
    # declare the vector for the updated angles
    T2_updated[zone] = Vector{Float64}()

    # get the reference node and index for the zone
    ref_node, shared_zone, shared_idx = Zonerefinfo[zone]
    ref_idx = nodename_nodeidx_map[zone][ref_node]

    # get the JuMP solution angle for the reference node
    current_ref_angle = value.(T2[zone][ref_idx])

    if shared_zone == nothing
      last_ref_angle = 0.0
    else
      # the shared_zone/shared_idx pair is the higher order shared zone and
      # node index for the reference node which is used to calculate the
      # adjustment needed in the current zone to make the angle values
      # the same
      last_ref_angle = T2_updated[shared_zone][shared_idx]
    end

    # calculate the difference or adjustment needed for each angle based
    # on the reference angle value in the higher zone order zone where the
    # reference node is shared and the current reference angle
    diff_angle = last_ref_angle - current_ref_angle
    println("zone $(zone), ref_node $(nodenames[zone][ref_idx]), last_ref_angle: $(last_ref_angle), current_ref_angle: $(current_ref_angle), diff_angle: $(diff_angle)")

    # setup for the next zone in the ordering by saving this reference angle
    # in order to calculate the next adjustment
    last_ref_angle = current_ref_angle

    # update every angle for the zone based on this adjustment factor
    for inode in 1:length(nodenames[zone])
      updated_angle = value.(T2[zone][inode]) + diff_angle
      append!(T2_updated[zone], updated_angle)
      println("zone $(zone), node $(nodenames[zone][inode]), original angle: $(value.(T2[zone][inode])), updated angle: $(T2_updated[zone][inode])")
    end
  end

  return T2_updated
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

function setup_estimate(measidxs, measidx_nodeidx_map, rmat, Ybus, Vnom)
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
  for inode = 1:nnode
    if inode in keys(Vnom)
      set_start_value.(v[inode], Vnom[inode][1])
    end
  end

  # Similar to magnitudes, for angles if we need the nominal values at all,
  # we should really use all of them to give us the best shot at finding the
  # correct solution for any model

  @variable(nlp,T[1:nnode])
  for inode = 1:nnode
    if inode in keys(Vnom)
      start = Vnom[inode][2]
      set_start_value.(T[inode], deg2rad(start))
      @NLconstraint(nlp, deg2rad(start-90.0) <= T[inode] <= deg2rad(start+90.0))
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


function perform_estimate(nlp, v, T)
  # call solver given everything is setup coming in
  @time optimize!(nlp)
  solution_summary(nlp, verbose=true)
  println("\nSolution v = $(value.(v))")
  println("\nSolution T = $(value.(T))")
end


function compare_estimate_magnitudes(v, nodenames, zone)
  inode = 0
  toterr = 0.0
  for row in CSV.File(string(test_dir, "/FPI_results_data.csv.", zone), header=false)
    inode += 1
    expected = row[2]
    solution = value.(v[inode])
    pererr = 100.0 * abs(solution - expected)/expected
    toterr += pererr
    println("$(nodenames[inode]) v exp: $(expected), sol: $(solution), %err: $(pererr)")
  end
  avgerr = toterr/inode
  println("*** Average v %err zone $(zone): $(avgerr)")
end


function compare_estimate_angles(T, nodenames, zone, Vnom)
  inode = 0
  toterr = 0.0
  for row in CSV.File(string(test_dir, "/FPI_results_data.csv.", zone), header=false)
    inode += 1
    expected = row[3]
    solution = T[inode]
    #solution = value.(T[inode])
    # hardwired logic for test case 180 degree Vnom angles
    if Vnom[inode][2] == 180
      solution -= pi
    end
    diff = abs(solution - expected)
    if diff < 1.0e-10
      pererr = 0.0
    else
      pererr = 100.0 * diff/expected
    end
    toterr += pererr
    println("$(nodenames[inode]) T exp: $(expected), sol: $(solution), %err: $(pererr)")
  end
  avgerr = toterr/inode
  println("*** Average T %err zone $(zone): $(avgerr)")
end


# Main

println("Start parsing input files...")

shared_nodenames = Dict()
for row in CSV.File(test_dir*"/sharednodes.csv", header=false)
  shared_nodenames[row[1]] = []
end

measidxs1 = Dict()
measidxs2 = Dict()
measidx1_nodeidx_map = Dict()
measidx2_nodeidx_map = Dict()
rmat1 = Dict()
rmat2 = Dict()
Ybus = Dict()
Ybusp = Dict()
Vnom = Dict()
nodenames = Dict()
nodename_nodeidx_map = Dict()
shared_nodeidx_measidx2_map = Dict()
measdata = Dict()

for zone = 0:5
  measidxs1[zone], measidxs2[zone], measidx1_nodeidx_map[zone], measidx2_nodeidx_map[zone], rmat1[zone], rmat2[zone], Ybus[zone], Ybusp[zone], Vnom[zone], nodenames[zone], nodename_nodeidx_map[zone], shared_nodeidx_measidx2_map[zone], measdata[zone] = get_input(zone, shared_nodenames)
end

Sharedmeas, SharedmeasAlt = setup_data_sharing(shared_nodenames, shared_nodeidx_measidx2_map)

# do the data structure initialization  for reference angle passing
Zoneorder, Zonerefinfo = setup_angle_passing(nodename_nodeidx_map, shared_nodenames)

println("Done parsing input files, start defining optimization problem...")

nlp1 = Dict()
nlp2 = Dict()
zvec1 = Dict()
zvec2 = Dict()
v1 = Dict()
v2 = Dict()
T1 = Dict()
T2 = Dict()

for zone = 0:5
  nlp1[zone], zvec1[zone], v1[zone], T1[zone] = setup_estimate(measidxs1[zone], measidx1_nodeidx_map[zone], rmat1[zone], Ybus[zone], Vnom[zone])
  nlp2[zone], zvec2[zone], v2[zone], T2[zone] = setup_estimate(measidxs2[zone], measidx2_nodeidx_map[zone], rmat2[zone], Ybus[zone], Vnom[zone])
end

println("\nDone defining optimization problem, start solving it...")

# assume all measurement_data files contain the same number of rows/timestamps
nrows = length(measdata[0])
println("number of timestamps to process: $(nrows)")

for row = 1:1 # first timestamp only
#for row = 1:nrows # all timestamps

  for zone = 0:5
    # This logic assumes that the order of measurement data (columns in a row)
    # is the same as measurement order (rows) in measurements.csv
    # If that's not the case, then this will require rework
    measurement = measdata[zone][row]
    for imeas in 1:length(measurement)-1
      set_value(zvec1[zone][imeas], measurement[imeas+1])
      set_value(zvec2[zone][imeas], measurement[imeas+1])
      println("*** zone $(zone) zvec1 and zvec2 initalization for measurement $(imeas): $(measurement[imeas+1])")
    end
    #println("*** Timestamp with measurements: $(measurement)\n")
  end

  # first optimization for each zone
  for zone = 0:5
    println("\n================================================================================")
    println("1st optimization for timestamp #$(row), zone: $(zone)\n")
    perform_estimate(nlp1[zone], v1[zone], T1[zone])
  end

  for zone = 0:5
    println("\n================================================================================")
    println("1st optimization magnitude comparison for timestamp #$(row), zone: $(zone)\n")
    compare_estimate_magnitudes(v1[zone], nodenames[zone], zone)
  end

  perform_data_sharing(Ybusp, Sharedmeas, SharedmeasAlt, measidxs2, v1, T1, zvec2)

  # second optimization after shared node data exchange
  for zone = 0:5
    println("\n================================================================================")
    println("2nd optimization for timestamp #$(row), zone: $(zone)\n")
    perform_estimate(nlp2[zone], v2[zone], T2[zone])
  end

  # perform reference angle passing to get the final angle results
  println("\n================================================================================")
  println("Reference angle passing:\n")
  T2_updated = perform_angle_passing(T2, Zoneorder, Zonerefinfo, nodename_nodeidx_map, nodenames)

  for zone = 0:5
    println("\n================================================================================")
    println("2nd optimization magnitude comparison for timestamp #$(row), zone: $(zone)\n")
    compare_estimate_magnitudes(v2[zone], nodenames[zone], zone)
  end

  for zone = 0:5
    println("\n================================================================================")
    println("2nd optimization angle comparison for timestamp #$(row), zone: $(zone)\n")
    compare_estimate_angles(T2_updated[zone], nodenames[zone], zone, Vnom[zone])
  end
end

