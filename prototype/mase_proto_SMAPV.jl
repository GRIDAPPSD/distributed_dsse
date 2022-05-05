#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays

test_dir = "mase_files_pu"


function goodbye()
  ccall(:jl_exit, Cvoid, (Int32,), 86)
end


function get_input(zone, shared_nodenames, Sharedalways_set)
  println("    Reading input files for zone: $(zone)")
  nodename_nodeidx_map = Dict()
  nodenames = Vector{String}()

  inode = 0
  for row in CSV.File(string(test_dir, "/nodelist.csv.", zone), header=false)
    inode += 1
    nodename_nodeidx_map[row[1]] = inode
    push!(nodenames, row[1])

    if haskey(shared_nodenames, row[1])
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
    if !haskey(measidxs1, stype)
      measidxs1[stype] = Vector{Int64}()
      measidxs2[stype] = Vector{Int64}()
    end
    imeas1 += 1
    append!(measidxs1[stype], imeas1)
    append!(measidxs2[stype], imeas1)

    measidx1_nodeidx_map[imeas1] = measidx2_nodeidx_map[imeas1] = nodename_nodeidx_map[row[3]]

    append!(rmat1, row[6]^2)
    append!(rmat2, row[6]^2)

    if haskey(shared_nodenames, row[3])
      nodeidx = nodename_nodeidx_map[row[3]]
      if !haskey(shared_nodeidx_measidx2_map, nodeidx)
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
      # create new measurement
      imeas2 += 1
      append!(measidxs2["vi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 0.01^2)
      println("***Adding missing vi measurement for zone: $(zone), shared node index: $(nodeidx), meas: $(imeas2)")
      push!(Sharedalways_set, (zone, imeas2))
    end
    if !PiFlag
      # create new measurement
      imeas2 += 1
      append!(measidxs2["Pi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 1.0^2)
      println("***Adding missing Pi measurement for zone: $(zone), shared node index: $(nodeidx), meas: $(imeas2)")
      push!(Sharedalways_set, (zone, imeas2))
    end
    if !QiFlag
      # create new measurement
      imeas2 += 1
      append!(measidxs2["Qi"], imeas2)
      measidx2_nodeidx_map[imeas2] = nodeidx
      append!(measidx_vec, imeas2)
      append!(rmat2, 0.5^2)
      println("***Adding missing Qi measurement zone: $(zone), for shared node index: $(nodeidx), meas: $(imeas2)")
      push!(Sharedalways_set, (zone, imeas2))
    end
  end

  println("    Total number of 1st estimate measurements: $(imeas1)")
  if haskey(measidxs1, "vi")
    println("    Number of vi measurements: $(length(measidxs1["vi"]))")
    println("    vi measurement indices: $(measidxs1["vi"])")
  end
  if haskey(measidxs1, "Ti")
    println("    Number of Ti measurements: $(length(measidxs1["Ti"]))")
    println("    Ti measurement indices: $(measidxs1["Ti"])")
  end
  if haskey(measidxs1, "Pi")
    println("    Number of Pi measurements: $(length(measidxs1["Pi"]))")
    println("    Pi measurement indices: $(measidxs1["Pi"])")
  end
  if haskey(measidxs1, "Qi")
    println("    Number of Qi measurements: $(length(measidxs1["Qi"]))")
    println("    Qi measurement indices: $(measidxs1["Qi"])")
  end

  println("    Total number of 2nd estimate measurements: $(imeas2)")
  if haskey(measidxs2, "vi")
    println("    Number of vi measurements: $(length(measidxs2["vi"]))")
    println("    vi measurement indices: $(measidxs2["vi"])")
  end
  if haskey(measidxs2, "Ti")
    println("    Number of Ti measurements: $(length(measidxs2["Ti"]))")
    println("    Ti measurement indices: $(measidxs2["Ti"])")
  end
  if haskey(measidxs2, "Pi")
    println("    Number of Pi measurements: $(length(measidxs2["Pi"]))")
    println("    Pi measurement indices: $(measidxs2["Pi"])")
  end
  if haskey(measidxs2, "Qi")
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
  for row in CSV.File(string(test_dir, "/vnom.csv.", zone))
    node = row[1]
    if haskey(nodename_nodeidx_map, node)
      Vnom[nodename_nodeidx_map[node]] = (row[2], row[3])
    end
  end
  println("    Vnom number of elements: $(length(Vnom))")

  measdata = CSV.File(string(test_dir, "/measurement_data.csv.", zone))

  return measidxs1, measidxs2, measidx1_nodeidx_map, measidx2_nodeidx_map, rmat1, rmat2, Ybus, Ybusp, Vnom, nodenames, nodename_nodeidx_map, shared_nodeidx_measidx2_map, measdata
end


function fake_variance_comparison()
  return true
end


function setup_data_sharing(shared_nodenames, shared_nodeidx_measidx2_map, Sharedalways_set)
  #Sharednodes = Dict()
  Sharedmeas = Dict()
  SharedmeasAlt = Dict()
  Secondestimate_set = Set()

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
      # determine whether data sharing should be done
      if (zone, imeas) in Sharedalways_set
        println("Always sharing data due to newly added measurement for destzone: $(zone), destmeas: $(imeas), sourcezone: $(value[2][1]), sourcenodeidx: $(value[2][2])")
      elseif fake_variance_comparison()
        println("Will be comparing variance, for now falling through to data sharing for destzone: $(zone), destmeas: $(imeas), sourcezone: $(value[2][1]), sourcenodeidx: $(value[2][2])")
      else
        # skip data sharing because it doesn't meet previous checks
        println("Skip sharing data for zone: $(zone), meas: $(imeas)")
        continue
      end

      Sharedmeas[zone][imeas] = value[2]
      if zone == 0
        SharedmeasAlt[zone][imeas] = value[2]
      else
        SharedmeasAlt[zone][imeas] = value[1]
      end
      push!(Secondestimate_set, zone)
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
      # determine whether data sharing should be done
      if (zone, imeas) in Sharedalways_set
        println("Always reverse sharing data due to newly added measurement for destzone: $(zone), destmeas: $(imeas), sourcezone: $(value[1][1]), sourcenodeidx: $(value[1][2])")
      elseif fake_variance_comparison()
        println("Will be comparing variance, for now falling through to reverse data sharing for destzone: $(zone), destmeas: $(imeas), sourcezone: $(value[1][1]), sourcenodeidx: $(value[1][2])")
      else
        # skip data sharing because it doesn't meet previous checks
        println("Skip sharing data for zone: $(zone), meas: $(imeas)")
        continue
      end

      Sharedmeas[zone][imeas] = value[1]
      if zone == 0
        SharedmeasAlt[zone][imeas] = value[1]
      else
        SharedmeasAlt[zone][imeas] = value[2]
      end
      push!(Secondestimate_set, zone)
    end
  end
  println("Shared nodenames dictionary: $(shared_nodenames)\n")
  #println("Sharednodes dictionary: $(Sharednodes)\n")
  println("Sharedmeas dictionary: $(Sharedmeas)\n")
  println("SharedmeasAlt dictionary: $(SharedmeasAlt)\n")
  println("Secondestimate_set: $(Secondestimate_set)\n")

  return Sharedmeas, SharedmeasAlt, Secondestimate_set
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

    V = Vector{ComplexF64}(undef, nnode)
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
      if haskey(measidxs2[zonedest], "vi") && measdest in measidxs2[zonedest]["vi"]
        println("vi measurement value sharing destzone: $(zonedest), destmeas: $(measdest), sourcezone: $(source[1]), sourcenode: $(source[2]), v1 sourcevalue: $(value.(v1[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], value.(v1[source[1]][source[2]]))
      elseif haskey(measidxs2[zonedest], "Ti") && measdest in measidxs2[zonedest]["Ti"]
        println("Ti measurement value sharing destzone: $(zonedest), destmeas: $(measdest), sourcezone: $(source[1]), sourcenode: $(source[2]), T1 sourcevalue: $(value.(T1[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], value.(T1[source[1]][source[2]]))
      elseif haskey(measidxs2[zonedest], "Pi") && measdest in measidxs2[zonedest]["Pi"]
        println("Pi measurement value sharing destzone: $(zonedest), destmeas: $(measdest), sourcezone: $(source[1]), sourcenode: $(source[2]), -S real sourcevalue: $(-1*real(S[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], -1*real(S[source[1]][source[2]]))
      elseif haskey(measidxs2[zonedest], "Qi") && measdest in measidxs2[zonedest]["Qi"]
        println("Qi measurement value sharing destzone: $(zonedest), destmeas: $(measdest), sourcezone: $(source[1]), sourcenode: $(source[2]), -S imag sourcevalue: $(-1*imag(S[source[1]][source[2]]))")
        set_value(zvec2[zonedest][measdest], -1*imag(S[source[1]][source[2]]))
      end
    end
  end
end


function close_to(value, check)
  return value-5.0<=check && value+5.0>=check
end


function get_phase(angle, phase_shift_flag)
  phase_shift = phase_shift_flag ? 30.0 : 0.0

  if close_to(angle, 0.0+phase_shift) || close_to(angle, 180.0+phase_shift)
    return "A"
  elseif close_to(angle, -120.0+phase_shift) || close_to(angle, 60.0+phase_shift)
    return "B"
  elseif close_to(angle, 120.0+phase_shift) || close_to(angle, -60.0+phase_shift)
    return "C"
  else
    println("WARNING: Vnom angle of $(angle) does not map to a phase, defaulting to phase A")
    return "A"
  end
end


function buildZonegraph(parzone, shared_nodenames, Zonenodes, Zonegraph)
  println("buildZonegraph parent zone: $(parzone)")
  for node in Zonenodes[parzone]
    for (zone, nodeidx) in shared_nodenames[node]
      if zone != parzone && !haskey(Zonegraph, zone)
        if !haskey(Zonegraph, parzone)
          Zonegraph[parzone] = Array{Tuple{Int64,Int64},1}()
        end

        buildZonegraph(zone, shared_nodenames, Zonenodes, Zonegraph)

        println("buildZonegraph Zonegraph[$(parzone)] add child zone: $(zone)")
        if haskey(Zonegraph, zone)
          push!(Zonegraph[parzone], (zone, length(Zonegraph[zone])))
        else
          push!(Zonegraph[parzone], (zone, 0))
        end
      end
    end
  end
end


function buildZoneorder(parzone, Zonegraph, Zoneorder)
  println("buildZoneorder parent zone: $(parzone)")
  if haskey(Zonegraph, parzone)
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


function setup_angle_passing(nodenames, nodename_nodeidx_map, Vnom, phase_set, phase_nodenames, phase_shared_nodenames)

  # first, determine if vnom angles are shifted 30 degrees based on the phase
  # A source node, thus allowing the phases of all other nodes to be determined
  phase_shift_flag = nothing

  for row in CSV.File(test_dir*"/sourcenodes.csv", header=false)
    node = row[1]
    for zone = 0:5
      if haskey(nodename_nodeidx_map[zone], node)
        angle = Vnom[zone][nodename_nodeidx_map[zone][node]][2]
        if close_to(angle, 0.0)
          phase_shift_flag = false
          break
        elseif close_to(angle, 30.0)
          phase_shift_flag = true
          break
        end
      end
    end
    if phase_shift_flag != nothing
      break
    end
  end

  if phase_shift_flag == nothing
    println("ERROR: could not determine the source node for phase A")
    goodbye()
  end

  # next, create phase_nodenames, phase_shared_nodenames, and phase_set
  # now that we know how to determine the phase
  for zone = 0:5
    phase_nodenames[zone] = Dict()

    for (nodeidx, (magnitude, angle)) in Vnom[zone]
      node = nodenames[zone][nodeidx]

      # determine the phase based on vnom angle
      phase = get_phase(angle, phase_shift_flag)
      if !haskey(phase_nodenames[zone], phase)
        phase_nodenames[zone][phase] = Vector{String}()
      end
      # add the node to the list of nodes for that phase
      push!(phase_nodenames[zone][phase], node)

      if haskey(shared_nodenames, node)
        if !haskey(phase_shared_nodenames, phase)
          phase_shared_nodenames[phase] = Dict()
        end
        phase_shared_nodenames[phase][node] = shared_nodenames[node]
      end

      push!(phase_set, phase)
    end
    #println("zone: $(zone), phase_nodenames: $(phase_nodenames[zone])")
  end
  #println("phase_set: $(phase_set)")
  #println("phase_shared_nodenames: $(phase_shared_nodenames)")

  # next, find and validate the system reference zone and node pair per
  # phase from source nodes
  sysref = Dict()
  for row in CSV.File(test_dir*"/sourcenodes.csv", header=false)
    node = row[1]
    for zone = 0:5
      if haskey(nodename_nodeidx_map[zone], node)
        # find the phase for the node
        for (phase, phasenodes) in phase_nodenames[zone]
          if node in phasenodes
            if haskey(sysref, phase)
              println("ERROR: found multiple source nodes for phase: $(phase)")
              goodbye()
            else
              sysref[phase] = (zone, node)
              break
            end
          end
        end
      end
    end
  end

  # check if there is a system reference for each phase that has nodes
  for phase in phase_set
    if !haskey(sysref, phase)
      println("ERROR: no system reference zone and node found for phase: $(phase)")
      goodbye()
    end
  end

  # check if the system reference for all phases are in the same zone
  checkzone = nothing
  for (phase, (refzone, refnode)) in sysref
    println("Found system reference for phase: $(phase), zone: $(refzone), node: $(refnode)")
    if checkzone == nothing
      checkzone = refzone
    elseif refzone != checkzone
      println("ERROR: system reference nodes not all in the same zone")
      goodbye()
    end
  end

  # next, determine the zone ordering
  Zoneorder = Dict()
  Zonerefinfo = Dict()

  for (phase, (refzone, refnode)) in sysref
    # build a graph of the zones linked to other zones by shared nodes
    # phase_shared_nodenames has the info needed to build this
    println("phase_shared_nodenames for phase $(phase): $(phase_shared_nodenames[phase])")

    # for each zone create a list of shared nodes
    Zonenodes = Dict()
    for (node, zonepairs) in phase_shared_nodenames[phase]
      for zonepair in zonepairs
        zone = zonepair[1]
        if !haskey(Zonenodes, zone)
          Zonenodes[zone] = Vector{String}()
        end
        push!(Zonenodes[zone], node)
      end
    end
    println("Shared nodes for phase $(phase) per zone, Zonenodes: $(Zonenodes)")

    Zonegraph = Dict()
    # invoke recursive function to build the graph of connected zones
    # pass the system reference zone and recursion will build the rest
    buildZonegraph(refzone, phase_shared_nodenames[phase], Zonenodes, Zonegraph)
    println("Connected zones graph, Zonegraph: $(Zonegraph)")

    Zoneorder[phase] = Vector{Int64}()
    append!(Zoneorder[phase], refzone) # system reference zone is always first
    # traverse the zone graph recursively starting from the system reference
    # zone to build the full zone ordering
    buildZoneorder(refzone, Zonegraph, Zoneorder[phase])
    println("Zone ordering vector for phase $(phase), Zoneorder: $(Zoneorder[phase])")

    # create a dictionary to quickly lookup order by zone
    iorder = 0
    ZoneorderDict = Dict()
    # iterate over Zoneorder backwards so the higher priority zones
    # get larger values
    for zone in Iterators.Reverse(Zoneorder[phase])
      iorder += 1
      ZoneorderDict[zone] = iorder
    end
    println("Zone ordering dictionary, ZoneorderDict: $(ZoneorderDict)")

    # finally, determine the zone reference nodes
    Zonerefinfo[phase] = Dict()

    for zone = 0:5
      if zone == refzone
        # for the system reference zone, the zone reference node is always
        # the system reference node
        Zonerefinfo[phase][zone] = (refnode, nothing, nothing)
      else
        # determine what the shared nodes are for this zone to find the one
        # that is the zone reference node
        # check each shared node to see which has the highest zone order for
        # the other zones where it is shared
        max_priority = 0
        max_node = max_zone = max_idx = nothing
        for node in Zonenodes[zone]
          # find other zone that node is shared with
          for (shared_zone, shared_idx) in phase_shared_nodenames[phase][node]
            if shared_zone!=zone && ZoneorderDict[shared_zone]>max_priority
              max_priority = ZoneorderDict[shared_zone]
              max_node = node
              max_zone = shared_zone
              max_idx = shared_idx
            end
          end
        end
        Zonerefinfo[phase][zone] = (max_node, max_zone, max_idx)
      end
    end
    println("Zone reference info for phase $(phase), Zonerefinfo: $(Zonerefinfo[phase])")
  end

  return Zoneorder, Zonerefinfo
end


function perform_angle_passing(T, Zoneorder, Zonerefinfo, nodename_nodeidx_map, phase_set, phase_nodenames)
  # store the updated angles after reference angle passing in a new data
  # structure because I can't update the JuMP T solution vector
  T_updated = Dict()
  for zone = 0:5
    # declare and allocate the per zone vectors for the updated angles
    T_updated[zone] = Vector{Float64}(undef, length(nodename_nodeidx_map[zone]))
  end

  for phase in phase_set
    for zone in Zoneorder[phase]
      # get the reference node and index for the zone along with the shared zone
      # and node index in that zone for the reference node
      ref_node, shared_zone, shared_idx = Zonerefinfo[phase][zone]
      ref_idx = nodename_nodeidx_map[zone][ref_node]

      # get the JuMP solution angle for the reference node
      current_ref_angle = value.(T[zone][ref_idx])

      if shared_zone == nothing
        last_ref_angle = 0.0
      else
        # the shared_zone/shared_idx pair is the higher order shared zone and
        # node index for the reference node which is used to calculate the
        # adjustment needed in the current zone to make the angle values
        # the same
        last_ref_angle = T_updated[shared_zone][shared_idx]
      end

      # calculate the difference or adjustment needed for each angle based
      # on the reference angle value in the higher zone order zone where the
      # reference node is shared and the current reference angle
      diff_angle = last_ref_angle - current_ref_angle
      println("phase $(phase), zone $(zone), ref_node $(nodenames[zone][ref_idx]), last_ref_angle: $(last_ref_angle), current_ref_angle: $(current_ref_angle), diff_angle: $(diff_angle)")

      # update every angle for nodes of the current phase in the zone based on
      # this adjustment factor
      for node in phase_nodenames[zone][phase]
        inode = nodename_nodeidx_map[zone][node]
        updated_angle = value.(T[zone][inode]) + diff_angle
        T_updated[zone][inode] = updated_angle
        println("phase $(phase), zone $(zone), node $(node), original angle: $(value.(T[zone][inode])), updated angle: $(T_updated[zone][inode])")
      end
    end
  end

  return T_updated
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
    if haskey(Vnom, inode)
      set_start_value.(v[inode], Vnom[inode][1])
    end
  end

  # Similar to magnitudes, for angles if we need the nominal values at all,
  # we should really use all of them to give us the best shot at finding the
  # correct solution for any model

  @variable(nlp,T[1:nnode])
  for inode = 1:nnode
    if haskey(Vnom, inode)
      start = Vnom[inode][2]
      set_start_value.(T[inode], deg2rad(start))
      @NLconstraint(nlp, deg2rad(start-90.0) <= T[inode] <= deg2rad(start+90.0))
    end
  end

  # Objective function formulation
  @NLexpression(nlp, vzi[i=1:nmeas], v[measidx_nodeidx_map[i]])

  # vi measurements
  if haskey(measidxs, "vi")
    @NLexpression(nlp, Visum, sum((zvec[i] - vzi[i])^2/rmat[i] for i in measidxs["vi"]))
  else
    @NLexpression(nlp, Visum, 0.0)
  end

  @NLexpression(nlp, Tzi[i=1:nmeas], T[measidx_nodeidx_map[i]])

  # Ti measurements
  if haskey(measidxs, "Ti")
    @NLexpression(nlp, Tisum, sum((zvec[i] - Tzi[i])^2/rmat[i] for i in measidxs["Ti"]))
  else
    @NLexpression(nlp, Tisum, 0.0)
  end

  # common terms for Pi and Qi measurements
  if haskey(measidxs, "Pi") || haskey(measidxs, "Qi")
    @NLexpression(nlp, Tzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], Tzi[i] - T[j])
    @NLexpression(nlp, Gzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], real.(Ybus[measidx_nodeidx_map[i]][j]))
    @NLexpression(nlp, Bzij[i=1:nmeas,j in keys(Ybus[measidx_nodeidx_map[i]])], imag.(Ybus[measidx_nodeidx_map[i]][j]))
  end

  # Pi measurements
  if haskey(measidxs, "Pi")
    @NLexpression(nlp, h_Pi[i in measidxs["Pi"]], vzi[i] * sum(v[j]*(Gzij[i,j]*cos(Tzij[i,j]) + Bzij[i,j]*sin(Tzij[i,j])) for j in keys(Ybus[measidx_nodeidx_map[i]])))
    @NLexpression(nlp, Pisum, sum((zvec[i] - h_Pi[i])^2/rmat[i] for i in measidxs["Pi"]))
  else
    @NLexpression(nlp, Pisum, 0.0)
  end

  # Qi measurements
  if haskey(measidxs, "Qi")
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


function compare_estimate_magnitudes(v, nodenames, FPIResults, zone, StatsMagnitude)
  toterr = 0.0
  nnode = length(v) # get number of nodes from # of v elements
  for inode = 1:nnode
    expected = FPIResults[inode][1]
    solution = value.(v[inode])
    diff = abs(solution - expected)

    if diff < StatsMagnitude[inode]["min"]
      StatsMagnitude[inode]["min"] = diff
    end
    if diff > StatsMagnitude[inode]["max"]
      StatsMagnitude[inode]["max"] = diff
    end
    #StatsMagnitude[inode]["sum"] += diff
    StatsMagnitude[inode]["sum"] = StatsMagnitude[inode]["sum"] + diff

    pererr = 100.0 * diff/expected
    toterr += pererr
    println("$(nodenames[inode]) v exp: $(expected), sol: $(solution), %err: $(pererr)")
  end
  avgerr = toterr/nnode
  println("*** Average v %err zone $(zone): $(avgerr)")
end


function compare_estimate_angles(T, nodenames, FPIResults, zone, Vnom, StatsAngle)
  toterr = 0.0
  span = 2*pi
  nnode = length(T) # get number of nodes from # of T elements
  for inode = 1:nnode
    expected = FPIResults[inode][2]
    solution = T[inode]
    #solution = value.(T[inode])
    # hardwired logic for test case 180 degree Vnom angles
    if Vnom[inode][2] == 180
      solution -= pi
    end

    diff = abs(solution - expected)
    if diff < StatsAngle[inode]["min"]
      StatsAngle[inode]["min"] = diff
    end
    if diff > StatsAngle[inode]["max"]
      StatsAngle[inode]["max"] = diff
    end
    #StatsAngle[inode]["sum"] += diff
    StatsAngle[inode]["sum"] = StatsAngle[inode]["sum"] + diff

    pererr = 100.0 * diff/span
    toterr += pererr
    println("$(nodenames[inode]) T exp: $(expected), sol: $(solution), 2*pi span %err: $(pererr)")
  end
  avgerr = toterr/nnode
  println("*** Average T 2*pi span %err zone $(zone): $(avgerr)")
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
Sharedalways_set = Set()

for zone = 0:5
  measidxs1[zone], measidxs2[zone], measidx1_nodeidx_map[zone], measidx2_nodeidx_map[zone], rmat1[zone], rmat2[zone], Ybus[zone], Ybusp[zone], Vnom[zone], nodenames[zone], nodename_nodeidx_map[zone], shared_nodeidx_measidx2_map[zone], measdata[zone] = get_input(zone, shared_nodenames, Sharedalways_set)
end

Sharedmeas, SharedmeasAlt, Secondestimate_set = setup_data_sharing(shared_nodenames, shared_nodeidx_measidx2_map, Sharedalways_set)

# do the data structure initialization for reference angle passing
phase_set = Set()
phase_nodenames = Dict()
phase_shared_nodenames = Dict()

Zoneorder, Zonerefinfo = setup_angle_passing(nodenames, nodename_nodeidx_map, Vnom, phase_set, phase_nodenames, phase_shared_nodenames)
#goodbye()

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

FPIResults = Dict()
StatsMagnitude1 = Dict()
StatsMagnitude2 = Dict()
StatsAngle1 = Dict()
StatsAngle2 = Dict()

for zone = 0:5
  FPIResults[zone] = Dict()
  StatsMagnitude1[zone] = Dict()
  StatsMagnitude2[zone] = Dict()
  StatsAngle1[zone] = Dict()
  StatsAngle2[zone] = Dict()
  nnode = length(Vnom[zone]) # get number of nodes from # of Vnom elements
  for row in CSV.File(string(test_dir, "/t_FPI_results_data.csv.", zone), header=true)
    timestamp = row[1]
    FPIResults[zone][timestamp] = Dict()
    for inode = 1:nnode
      FPIResults[zone][timestamp][inode] = (row[2*inode], row[2*inode+1])
    end
  end
  for inode = 1:nnode
    StatsMagnitude1[zone][inode] = Dict()
    StatsMagnitude1[zone][inode]["min"] = typemax(Float64)
    StatsMagnitude1[zone][inode]["max"] = 0.0
    StatsMagnitude1[zone][inode]["sum"] = 0.0
    StatsMagnitude2[zone][inode] = Dict()
    StatsMagnitude2[zone][inode]["min"] = typemax(Float64)
    StatsMagnitude2[zone][inode]["max"] = 0.0
    StatsMagnitude2[zone][inode]["sum"] = 0.0
    StatsAngle1[zone][inode] = Dict()
    StatsAngle1[zone][inode]["min"] = typemax(Float64)
    StatsAngle1[zone][inode]["max"] = 0.0
    StatsAngle1[zone][inode]["sum"] = 0.0
    StatsAngle2[zone][inode] = Dict()
    StatsAngle2[zone][inode]["min"] = typemax(Float64)
    StatsAngle2[zone][inode]["max"] = 0.0
    StatsAngle2[zone][inode]["sum"] = 0.0
  end
end

println("\nDone defining optimization problem, start solving it...")

# assume all measurement_data files contain the same number of rows/timestamps
nrows = length(measdata[0])
println("number of timestamps to process: $(nrows)")

ntimestamps = 0
for row = 1:1 # first timestamp only
#for row = 1:nrows # all timestamps
  global ntimestamps += 1

  for zone = 0:5
    # This logic assumes that the order of measurement data (columns in a row)
    # is the same as measurement order (rows) in measurements.csv
    # If that's not the case, then this will require rework
    measurement = measdata[zone][row]
    for imeas in 1:length(measurement)-1
      set_value(zvec1[zone][imeas], measurement[imeas+1])
      set_value(zvec2[zone][imeas], measurement[imeas+1])
      #println("*** zone $(zone) zvec1 and zvec2 initalization for measurement $(imeas): $(measurement[imeas+1])")
    end
    #println("*** Timestamp with measurements: $(measurement)\n")
  end

  # first optimization for each zone
  for zone = 0:5
    println("\n================================================================================")
    println("1st optimization for timestamp #$(row), zone: $(zone)\n")
    perform_estimate(nlp1[zone], v1[zone], T1[zone])
  end

  # perform reference angle passing to get the final angle results
  println("\n================================================================================")
  println("Reference angle passing:")
  T1_updated = perform_angle_passing(T1, Zoneorder, Zonerefinfo, nodename_nodeidx_map, phase_set, phase_nodenames)

  for zone = 0:5
    timestamp = measdata[zone][row][1]
    println("\n================================================================================")
    println("1st optimization magnitude comparison for timestamp $(timestamp), zone: $(zone)\n")
    compare_estimate_magnitudes(v1[zone], nodenames[zone], FPIResults[zone][timestamp], zone, StatsMagnitude1[zone])
  end

  for zone = 0:5
    timestamp = measdata[zone][row][1]
    println("\n================================================================================")
    println("1st optimization angle comparison for timestamp $(timestamp), zone: $(zone)\n")
    compare_estimate_angles(T1_updated[zone], nodenames[zone], FPIResults[zone][timestamp], zone, Vnom[zone], StatsAngle1[zone])
  end

  perform_data_sharing(Ybusp, Sharedmeas, SharedmeasAlt, measidxs2, v1, T1, zvec2)

  # second optimization after shared node data exchange
  for zone in Secondestimate_set
    println("\n================================================================================")
    println("2nd optimization for timestamp #$(row), zone: $(zone)\n")
    perform_estimate(nlp2[zone], v2[zone], T2[zone])
  end

  # perform reference angle passing to get the final angle results
  println("\n================================================================================")
  println("Reference angle passing:")
  T2_updated = perform_angle_passing(T2, Zoneorder, Zonerefinfo, nodename_nodeidx_map, phase_set, phase_nodenames)

  for zone = 0:5
    timestamp = measdata[zone][row][1]
    println("\n================================================================================")
    println("2nd optimization magnitude comparison for timestamp $(timestamp), zone: $(zone)\n")
    compare_estimate_magnitudes(v2[zone], nodenames[zone], FPIResults[zone][timestamp], zone, StatsMagnitude2[zone])
  end

  for zone = 0:5
    timestamp = measdata[zone][row][1]
    println("\n================================================================================")
    println("2nd optimization angle comparison for timestamp $(timestamp), zone: $(zone)\n")
    compare_estimate_angles(T2_updated[zone], nodenames[zone], FPIResults[zone][timestamp], zone, Vnom[zone], StatsAngle2[zone])
  end
end

mag_max_max = 0.0
mag_max_max_zone = ""
mag_max_max_node = ""
mag_max_mean = 0.0
mag_max_mean_zone = ""
mag_max_mean_node = ""
for zone = 0:5
  println("\n2nd optimization magnitude min, max, mean difference stats for zone: $(zone), # of timestamps: $(ntimestamps):")
  for inode = 1:length(nodenames[zone])
    global mag_max_max, mag_max_max_zone, mag_max_max_node
    global mag_max_mean, mag_max_mean_zone, mag_max_mean_node
    mean = StatsMagnitude2[zone][inode]["sum"]/ntimestamps
    println("""  $(nodenames[zone][inode]): $(StatsMagnitude2[zone][inode]["min"]), $(StatsMagnitude2[zone][inode]["max"]), $(mean)""")
    if StatsMagnitude2[zone][inode]["max"] > mag_max_max
      mag_max_max = StatsMagnitude2[zone][inode]["max"]
      mag_max_max_zone = zone
      mag_max_max_node = nodenames[zone][inode]
    end
    if mean > mag_max_mean
      mag_max_mean = mean
      mag_max_mean_zone = zone
      mag_max_mean_node = nodenames[zone][inode]
    end
  end
end

angle_max_max = 0.0
angle_max_max_zone = ""
angle_max_max_node = ""
angle_max_mean = 0.0
angle_max_mean_zone = ""
angle_max_mean_node = ""
for zone = 0:5
  println("\n2nd optimization angle min, max, mean difference stats for zone: $(zone), # of timestamps: $(ntimestamps):")
  for inode = 1:length(nodenames[zone])
    global angle_max_max, angle_max_max_zone, angle_max_max_node
    global angle_max_mean, angle_max_mean_zone, angle_max_mean_node
    mean = StatsAngle2[zone][inode]["sum"]/ntimestamps
    println("""  $(nodenames[zone][inode]): $(StatsAngle2[zone][inode]["min"]), $(StatsAngle2[zone][inode]["max"]), $(mean)""")
    if StatsAngle2[zone][inode]["max"] > angle_max_max
      angle_max_max = StatsAngle2[zone][inode]["max"]
      angle_max_max_zone = zone
      angle_max_max_node = nodenames[zone][inode]
    end
    if mean > angle_max_mean
      angle_max_mean = mean
      angle_max_mean_zone = zone
      angle_max_mean_node = nodenames[zone][inode]
    end
  end
end

# convert to degrees because it's easier when dealing with small values
angle_max_max = rad2deg(angle_max_max)
angle_max_mean = rad2deg(angle_max_mean)

println("\n2nd optimization magnitude max difference zone: $(mag_max_max_zone), node: $(mag_max_max_node), value: $(mag_max_max)")
println("2nd optimization magnitude max mean difference zone: $(mag_max_mean_zone), node: $(mag_max_mean_node), value: $(mag_max_mean)")
println("2nd optimization angle max difference zone: $(angle_max_max_zone), node: $(angle_max_max_node), value: $(angle_max_max)")
println("2nd optimization angle max mean difference zone: $(angle_max_mean_zone), node: $(angle_max_mean_node), value: $(angle_max_mean)")

