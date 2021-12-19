#!/usr/bin/env julia 
using JuMP
using Ipopt
using CSV
using SparseArrays


function parse_measurements()
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

  return nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs, Pi_measidxs, Qi_measidxs, YbusG, YbusB
end


# Main

#nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10,"max_iter"=>80))
nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-1, "max_iter"=>10000))

println("Start parsing files...")

nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs, Pi_measidxs, Qi_measidxs, YbusG, YbusB = parse_measurements()

println("Done parsing files, start defining optimization problem...")

# declare z here then set values later when iterating over
# measurement_data.csv
@NLparameter(nlp,zvec[i=1:nmeas] == 0)

v_target_solution = [2280.1462224653337216, 2370.9304940663460002, 2291.1990472686602516, 2279.9611038436473791, 2370.9035460218810840, 2291.0555032704160112, 2362.1820781162941785, 2337.7110672203370996, 2362.1808119583511143, 2337.7095243841390584, 2277.5348581324765291, 2290.0280729393575712, 2334.7435197623253771, 2373.6298734831798356, 2340.0148421246135513, 2316.8779210656639407, 2371.6005614140126454, 2322.1821129960808321, 2269.9339329437420929, 2401.5701366959156076, 2401.6286457108431023, 2401.5788837522886752, 2270.2660272739585707, 2372.2073571348391852, 2286.7886768931412007, 2288.8924150224479490, 2330.3982535612799438, 2371.3599533499755125, 2337.0292728757217446, 2280.1464074987165986, 2370.9306694021715884, 2291.1992213259291020, 66395.3000000000029104, 66395.2999999999883585, 66395.2999999999883585, 2401.7203364564375079, 2401.7385918339664386, 2401.7263986362404466, 265.3904791881118967, 270.9028417125392707, 266.9026023653719903]
T_target_solution = [-0.0480446300322439, -2.1104985943393637, 2.0545258974007106, -0.0480650512861194, -2.1104993103261331, 2.0545086644810366, -2.1086030027029521, 2.0722242874676491, -2.1086034364935085, 2.0722245481961461, -0.0480391337561771, 2.0529656575179640, -0.0229219642187778, -2.1069033570263000, 2.0719690478578934, -0.0311579982700007, -2.1080452706239958, 2.0654498954911631, -0.0476869700719445, -0.0000394796490696, -2.0944289185520728, 2.0943476493206190, -0.0494099448676698, -2.1112820332230458, 2.0554359284967827, 2.0509730331871685, -0.0234745221225297, -2.1072314747030054, 2.0717780313744405, -0.0480446535613301, -2.1104986233690131, 2.0545258706446399, 0.5235987755982988, -1.5707963267948963, 2.6179938779914944, -0.0000273674186166, -2.0944166624406009, 2.0943655287652119, -0.0297801199764095, -2.1114943921009415, 2.0673840447661922]

#@variable(nlp,-0.5 <= v[1:nnode] <= 1.5,start=1.0)
#@variable(nlp,-pi <= T[1:nnode] <= pi,start=0.0)
#@variable(nlp,v[1:nnode],start=1.0)
#@variable(nlp,v[1:nnode],start=2401.6)
#@variable(nlp,T[1:nnode],start=0.0)
@variable(nlp,v[1:nnode])
set_start_value.(v, v_target_solution)
@variable(nlp,T[1:nnode])
set_start_value.(T, T_target_solution)

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

#v_target_solution = [2280.1462224653337216, 2370.9304940663460002, 2291.1990472686602516, 2279.9611038436473791, 2370.9035460218810840, 2291.0555032704160112, 2362.1820781162941785, 2337.7110672203370996, 2362.1808119583511143, 2337.7095243841390584, 2277.5348581324765291, 2290.0280729393575712, 2334.7435197623253771, 2373.6298734831798356, 2340.0148421246135513, 2316.8779210656639407, 2371.6005614140126454, 2322.1821129960808321, 2269.9339329437420929, 2401.5701366959156076, 2401.6286457108431023, 2401.5788837522886752, 2270.2660272739585707, 2372.2073571348391852, 2286.7886768931412007, 2288.8924150224479490, 2330.3982535612799438, 2371.3599533499755125, 2337.0292728757217446, 2280.1464074987165986, 2370.9306694021715884, 2291.1992213259291020, 66395.3000000000029104, 66395.2999999999883585, 66395.2999999999883585, 2401.7203364564375079, 2401.7385918339664386, 2401.7263986362404466, 265.3904791881118967, 270.9028417125392707, 266.9026023653719903]
#T_target_solution = [-0.0480446300322439, -2.1104985943393637, 2.0545258974007106, -0.0480650512861194, -2.1104993103261331, 2.0545086644810366, -2.1086030027029521, 2.0722242874676491, -2.1086034364935085, 2.0722245481961461, -0.0480391337561771, 2.0529656575179640, -0.0229219642187778, -2.1069033570263000, 2.0719690478578934, -0.0311579982700007, -2.1080452706239958, 2.0654498954911631, -0.0476869700719445, -0.0000394796490696, -2.0944289185520728, 2.0943476493206190, -0.0494099448676698, -2.1112820332230458, 2.0554359284967827, 2.0509730331871685, -0.0234745221225297, -2.1072314747030054, 2.0717780313744405, -0.0480446535613301, -2.1104986233690131, 2.0545258706446399, 0.5235987755982988, -1.5707963267948963, 2.6179938779914944, -0.0000273674186166, -2.0944166624406009, 2.0943655287652119, -0.0297801199764095, -2.1114943921009415, 2.0673840447661922]

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


  v_diff_solution = value.(v) - v_target_solution
  T_diff_solution = value.(T) - T_target_solution

  println("\nv target solution difference:  $(v_diff_solution)")
  println("\nv maximum solution difference:  $(findmax(v_diff_solution))")

  println("\nT target solution difference:  $(T_diff_solution)")
  println("\nT maximum solution difference:  $(findmax(T_diff_solution))")
end

