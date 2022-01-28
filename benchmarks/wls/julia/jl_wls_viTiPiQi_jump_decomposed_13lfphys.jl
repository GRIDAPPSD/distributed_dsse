#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays


function parse_measurements()
  nodename_nodeidx_map = Dict()
  nodeidx_nodename_map = Dict()
  nnode = 0
  for row in CSV.File("test/nodelist.csv", header=false)
    nnode += 1
    nodename_nodeidx_map[row[1]] = nnode
    nodeidx_nodename_map[nnode] = row[1]
  end
  println("    Total number of nodes: $(nnode)")

  # figure out how many total measurements we have for allocating arrays
  nmeas = (open("test/measurements.csv") |> readlines |> length) -1
  println("    Total number of measurements: $(nmeas)")

  vi_measidxs = Vector{Int64}()
  Ti_measidxs = Vector{Int64}()
  Pi_measidxs = Vector{Int64}()
  Qi_measidxs = Vector{Int64}()
  measidx_nodeidx_map = Dict()
  rmat = Array{Float64}(undef, nmeas)

  imeas = 0
  for row in CSV.File("test/measurements.csv")
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
    elseif row[1] == "Pi"
      supported = true
      imeas += 1
      append!(Pi_measidxs, imeas)
    elseif row[1] == "Qi"
      supported = true
      imeas += 1
      append!(Qi_measidxs, imeas)
    end

    if supported
      rmat[imeas] = row[6]^2
      measidx_nodeidx_map[imeas] = nodename_nodeidx_map[row[3]]
    end
  end
  println("    Number of vi measurements: $(length(vi_measidxs))")
  println("    Number of Ti measurements: $(length(Ti_measidxs))")
  println("    Number of Pi measurements: $(length(Pi_measidxs))")
  println("    Number of Qi measurements: $(length(Qi_measidxs))")
  println("    vi measurement indices: $(vi_measidxs)")
  println("    Ti measurement indices: $(Ti_measidxs)")
  println("    Pi measurement indices: $(Pi_measidxs)")
  println("    Qi measurement indices: $(Qi_measidxs)")
  println("    measurement node index map: $(measidx_nodeidx_map)")

  #Ybus = Dict()
  YbusG = Dict()
  YbusB = Dict()
  for row in CSV.File("test/ysparse.csv")
    if !haskey(YbusG, row[1])
      #Ybus[row[1]] = Dict()
      YbusG[row[1]] = Dict()
      YbusB[row[1]] = Dict()
    end
    # must construct full Ybus, not just lower diagonal elements
    #Ybus[row[1]][row[2]] = Ybus[row[2]][row[1]] = complex(row[3], row[4])
    YbusG[row[1]][row[2]] = YbusG[row[2]][row[1]] = row[3]
    YbusB[row[1]][row[2]] = YbusB[row[2]][row[1]] = row[4]
  end

  return nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs, Pi_measidxs, Qi_measidxs, YbusG, YbusB, nodeidx_nodename_map
end


# Main

#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"max_iter"=>80))
#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-1, "max_iter"=>10000))
#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-1,"max_iter"=>10000))
#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-1,"max_iter"=>10000))
#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-9,"max_iter"=>10000))
nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"acceptable_tol"=>1e-10))
#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-14,"acceptable_tol"=>1e-14,"max_iter"=>100000))

println("Start parsing files...")

nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs, Pi_measidxs, Qi_measidxs, YbusG, YbusB, nodeidx_nodename_map = parse_measurements()

println("Done parsing files, start defining optimization problem...")

# declare z here then set values later when iterating over
# measurement_data.csv
@NLparameter(nlp,zvec[i=1:nmeas] == 0)

v_target_solution = [2283.0834341589052201,2372.3754031263556499,2294.3624201193733825,2282.8985641419048989,2372.3484716871507771,2294.2190783191026640,2363.2752997397833497,2339.3809287173489793,2361.2315564450464080,2336.8924815796344774,2280.4987463440056672,2293.2337589443063735,2336.1909734539181045,2374.6161345817708934,2341.6688920407827936,2318.8179479999921568,2372.7592114495701026,2324.3778369779993227,2273.0063120857521426,2401.4131985750614149,2401.5134749777644174,2401.4247885560480427,2273.3742005452431840,2373.6699522389599224,2290.0518874935391977,2292.1405804143196292,2331.8752568378727119,2372.3604465902581069,2338.7013584566648206,2283.0836188232692621,2372.3755788026765003,2294.3625947736331909,66395.3000000000174623,66395.2999999999883585,66395.2999999999883585,2401.7203456675565576,2401.7385699396149903,2401.7264235380616810,265.5451524017101974,271.0012011180325544,267.0795161521803038]
T_target_solution = [-0.0487008581519078,-2.1110176173411861,2.0538786105973394,-0.0487212231106527,-2.1110183324034084,2.0538614270470061,-2.1089698863920869,2.0717849099347547,-2.1096690990054103,2.0722047197955713,-0.0487015503792061,2.0523279899930449,-0.0233251756205417,-2.1072557600226691,2.0715461941386337,-0.0316480734366957,-2.1084606908889381,2.0649449934416371,-0.0483711665236682,-0.0000521372337437,-2.0944418330827110,2.0943289314156139,-0.0500893376253135,-2.1118074287584010,2.0547562138638016,2.0503481067922014,-0.0239046059606650,-2.1076067335546274,2.0713291512159540,-0.0487008796533500,-2.1110176441267274,2.0538785857932265,0.5235987755982987,-1.5707963267948968,2.6179938779914944,-0.0000273642134291,-2.0944166946151137,2.0943655274575939,-0.0303208337814400,-2.1119850227417052,2.0668226042254458]

#@variable(nlp,v[1:nnode],start=66395.3)
@variable(nlp,v[1:nnode],start=2401.6)
#@variable(nlp,v[1:nnode])
#for i = 1:nnode
#  if i>=33 && i<=35
#    set_start_value.(v[i], 66395.3)
#  elseif i>=39 && i<=41
#    set_start_value.(v[i], 277.13)
#  else
#    set_start_value.(v[i], 2401.6)
#  end
#end

#@variable(nlp,T[1:nnode],start=0.0)
#@variable(nlp,-pi <= T[1:nnode] <= pi)
@variable(nlp,T[1:nnode])
for i = 1:nnode
  node = nodeidx_nodename_map[i]
  if endswith(node, ".1")
    set_start_value.(T[i], 0.0)
    #set_start_value.(T[i], deg2rad(30.0))
    #@NLconstraint(nlp, T[i] >= deg2rad(30.0-90.0))
    #@NLconstraint(nlp, T[i] <= deg2rad(30.0+90.0))
  elseif endswith(node, ".2")
    set_start_value.(T[i], deg2rad(-120.0))
    #set_start_value.(T[i], deg2rad(-90.0))
    #@NLconstraint(nlp, T[i] >= deg2rad(-90.0-90.0))
    #@NLconstraint(nlp, T[i] <= deg2rad(-90.0+90.0))
  else
    set_start_value.(T[i], deg2rad(120.0))
    #set_start_value.(T[i], deg2rad(150.0))
    #@NLconstraint(nlp, T[i] >= deg2rad(150.0-90.0))
    #@NLconstraint(nlp, T[i] <= deg2rad(150.0+90.0))
  end
end

set_start_value.(v[33], 66395.3)
set_start_value.(v[34], 66395.3)
set_start_value.(v[35], 66395.3)
@NLconstraint(nlp, [i=33:35], v[i] == 66395.3)

set_start_value.(T[33], deg2rad(30.0))
set_start_value.(T[34], deg2rad(-90.0))
set_start_value.(T[35], deg2rad(150.0))
@NLconstraint(nlp, T[33] == deg2rad(30.0))
@NLconstraint(nlp, T[34] == deg2rad(-90.0))
@NLconstraint(nlp, T[35] == deg2rad(150.0))

#@variable(nlp,v[1:nnode])
#set_start_value.(v, v_target_solution)
#@variable(nlp,T[1:nnode])
#set_start_value.(T, T_target_solution)

@NLexpression(nlp, vzi[i=1:nmeas], v[measidx_nodeidx_map[i]])

@NLexpression(nlp, Visum, sum((zvec[i] - vzi[i])^2/rmat[i] for i in vi_measidxs))

@NLexpression(nlp, Tzi[i=1:nmeas], T[measidx_nodeidx_map[i]])

@NLexpression(nlp, Tisum, sum((zvec[i] - Tzi[i])^2/rmat[i] for i in Ti_measidxs))

@NLexpression(nlp, Tzij[i=1:nmeas,j in keys(YbusG[measidx_nodeidx_map[i]])], Tzi[i] - T[j])
@NLexpression(nlp, Gzij[i=1:nmeas,j in keys(YbusG[measidx_nodeidx_map[i]])], YbusG[measidx_nodeidx_map[i]][j])
@NLexpression(nlp, Bzij[i=1:nmeas,j in keys(YbusB[measidx_nodeidx_map[i]])], YbusB[measidx_nodeidx_map[i]][j])

@NLexpression(nlp, h_Pi[i in Pi_measidxs], vzi[i] * sum(v[j]*(Gzij[i,j]*cos(Tzij[i,j]) + Bzij[i,j]*sin(Tzij[i,j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))) 
@NLexpression(nlp, Pisum, sum((zvec[i] - h_Pi[i])^2/rmat[i] for i in Pi_measidxs))

@NLexpression(nlp, h_Qi[i in Qi_measidxs], vzi[i] * sum(v[j]*(Gzij[i,j]*sin(Tzij[i,j]) - Bzij[i,j]*cos(Tzij[i,j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))) 
@NLexpression(nlp, Qisum, sum((zvec[i] - h_Qi[i])^2/rmat[i] for i in Qi_measidxs))

@NLobjective(nlp, Min, Visum + Tisum + Pisum + Qisum)

#@NLobjective(nlp, Min,
#  sum((zvec[i] - v[measidx_nodeidx_map[i]])^2/rmat[i] for i in vi_measidxs) +
#
#  sum((zvec[i] - T[measidx_nodeidx_map[i]])^2/rmat[i] for i in Ti_measidxs) +
#
#  sum((zvec[i] - (v[measidx_nodeidx_map[i]] * sum(v[j]*(YbusG[measidx_nodeidx_map[i]][j]*cos(T[measidx_nodeidx_map[i]] - T[j]) + YbusB[measidx_nodeidx_map[i]][j]*sin(T[measidx_nodeidx_map[i]] - T[j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Pi_measidxs) +
#
#  sum((zvec[i] - (v[measidx_nodeidx_map[i]] * sum(v[j]*(YbusG[measidx_nodeidx_map[i]][j]*sin(T[measidx_nodeidx_map[i]] - T[j]) - YbusB[measidx_nodeidx_map[i]][j]*cos(T[measidx_nodeidx_map[i]] - T[j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Qi_measidxs))


println("Done with defining optimization problem, start solving it...")

# process each timestamp getting measurement data and calling solver

for row in CSV.File("test/measurement_data.csv")
  println("\n================================================================================\n")
  # This logic assumes that the order of measurements (columns in a row)
  # in meaurement_data.csv is the same as measurement order (rows)
  # in measurements.csv
  # If that's not the case, then this will require some fanciness with
  # dictionaries to assign zvec values to the appropriate index
  for imeas in 1:nmeas
    set_value(zvec[imeas], row[imeas+1])
  end
  println("*** Timestamp with measurements: $(row)\n")

  @time optimize!(nlp)
  solution_summary(nlp, verbose=true)
  println("v = $(value.(v))")
  println("T = $(value.(T))")


  v_diff_solution = abs.(value.(v) - v_target_solution)
  T_diff_solution = abs.(value.(T) - T_target_solution)

  println("\nv target solution difference:  $(v_diff_solution)")
  println("\nv maximum solution difference:  $(findmax(v_diff_solution))")

  println("\nT target solution difference:  $(T_diff_solution)")
  println("\nT maximum solution difference:  $(findmax(T_diff_solution))")

  println("\nAndy special: expected target vs. solution side-by-side and diff...")
  for i = 1:nnode
    println("$(nodeidx_nodename_map[i]) v exp: $(v_target_solution[i]), sol: $(value.(v[i])), diff: $(abs(v_target_solution[i]-value.(v[i])))")
  end
  println("")
  for i = 1:nnode
    println("$(nodeidx_nodename_map[i]) T exp: $(T_target_solution[i]), sol: $(value.(T[i])), diff: $(abs(T_target_solution[i]-value.(T[i])))")
  end
end

