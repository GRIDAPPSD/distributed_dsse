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

# For some reason I'm not aware of yet at least, this code is way slower
# when it's wrapped by a function even though my understanding is that
# Julia compiles functions so they are supposed to perform better than
# code in the main/global scope. So keep this in main until I figure it out.
#function optimize()
  println("Start parsing files...")

  nnode, nz, rmat, zidx_nodeidx_map, vi_zidxs, Ti_zidxs, Pi_zidxs, Qi_zidxs, Ybus = parse_measurements()

  println("Done parsing files, start defining optimization problem...")

  # just hardwire initial x starting point, x0, for now
  x0 = Array{Float64}(undef, 2*nnode)
  for i = 1:nnode
  #  x0[i] = 1.0 # initial magnitudes
    x0[i] = 2401.6 # initial magnitudes
  end
  for i = nnode+1:2*nnode
    x0[i] = 0.0 # initial angles
  end
  x0 = [2283.0834341589052201,2372.3754031263556499,2294.3624201193733825,2282.8985641419048989,2372.3484716871507771,2294.2190783191026640,2363.2752997397833497,2339.3809287173489793,2361.2315564450464080,2336.8924815796344774,2280.4987463440056672,2293.2337589443063735,2336.1909734539181045,2374.6161345817708934,2341.6688920407827936,2318.8179479999921568,2372.7592114495701026,2324.3778369779993227,2273.0063120857521426,2401.4131985750614149,2401.5134749777644174,2401.4247885560480427,2273.3742005452431840,2373.6699522389599224,2290.0518874935391977,2292.1405804143196292,2331.8752568378727119,2372.3604465902581069,2338.7013584566648206,2283.0836188232692621,2372.3755788026765003,2294.3625947736331909,66395.3000000000174623,66395.2999999999883585,66395.2999999999883585,2401.7203456675565576,2401.7385699396149903,2401.7264235380616810,265.5451524017101974,271.0012011180325544,267.0795161521803038,-0.0487008581519078,-2.1110176173411861,2.0538786105973394,-0.0487212231106527,-2.1110183324034084,2.0538614270470061,-2.1089698863920869,2.0717849099347547,-2.1096690990054103,2.0722047197955713,-0.0487015503792061,2.0523279899930449,-0.0233251756205417,-2.1072557600226691,2.0715461941386337,-0.0316480734366957,-2.1084606908889381,2.0649449934416371,-0.0483711665236682,-0.0000521372337437,-2.0944418330827110,2.0943289314156139,-0.0500893376253135,-2.1118074287584010,2.0547562138638016,2.0503481067922014,-0.0239046059606650,-2.1076067335546274,2.0713291512159540,-0.0487008796533500,-2.1110176441267274,2.0538785857932265,0.5235987755982987,-1.5707963267948968,2.6179938779914944,-0.0000273642134291,-2.0944166946151137,2.0943655274575939,-0.0303208337814400,-2.1119850227417052,2.0668226042254458]

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
  #Vi(x) = sum((zvec[zidx] - x[zidx_nodeidx_map[zidx]])^2/rmat[zidx] for zidx in vi_zidxs)
  #Ti(x) = sum((zvec[zidx] - x[zidx_nodeidx_map[zidx]+nnode])^2/rmat[zidx] for zidx in Ti_zidxs)


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
  #target_solution = [1.0000, 0.971576218875089, 0.956643047805180, 0.950535446461549, 0.0, 0.0, 0.0, 0.0] # FPI solver 10,000 iterations
  # MATLAB FPI solver solution for the 3p6 test case 
  #target_solution = [1.0, 1.0, 1.0, 0.9712111193225086, 0.9410836694466717, 0.9223652753951702, 0.0, -2.0943951023931957, 2.0943951023931953, -0.0375155389352858, -2.1362637169183625, 2.0491375942741321] 
  # MATLAB FPI solver solution for the 13assets_lf physical units test case 
  target_solution = [2283.0834341589052201,2372.3754031263556499,2294.3624201193733825,2282.8985641419048989,2372.3484716871507771,2294.2190783191026640,2363.2752997397833497,2339.3809287173489793,2361.2315564450464080,2336.8924815796344774,2280.4987463440056672,2293.2337589443063735,2336.1909734539181045,2374.6161345817708934,2341.6688920407827936,2318.8179479999921568,2372.7592114495701026,2324.3778369779993227,2273.0063120857521426,2401.4131985750614149,2401.5134749777644174,2401.4247885560480427,2273.3742005452431840,2373.6699522389599224,2290.0518874935391977,2292.1405804143196292,2331.8752568378727119,2372.3604465902581069,2338.7013584566648206,2283.0836188232692621,2372.3755788026765003,2294.3625947736331909,66395.3000000000174623,66395.2999999999883585,66395.2999999999883585,2401.7203456675565576,2401.7385699396149903,2401.7264235380616810,265.5451524017101974,271.0012011180325544,267.0795161521803038,-0.0487008581519078,-2.1110176173411861,2.0538786105973394,-0.0487212231106527,-2.1110183324034084,2.0538614270470061,-2.1089698863920869,2.0717849099347547,-2.1096690990054103,2.0722047197955713,-0.0487015503792061,2.0523279899930449,-0.0233251756205417,-2.1072557600226691,2.0715461941386337,-0.0316480734366957,-2.1084606908889381,2.0649449934416371,-0.0483711665236682,-0.0000521372337437,-2.0944418330827110,2.0943289314156139,-0.0500893376253135,-2.1118074287584010,2.0547562138638016,2.0503481067922014,-0.0239046059606650,-2.1076067335546274,2.0713291512159540,-0.0487008796533500,-2.1110176441267274,2.0538785857932265,0.5235987755982987,-1.5707963267948968,2.6179938779914944,-0.0000273642134291,-2.0944166946151137,2.0943655274575939,-0.0303208337814400,-2.1119850227417052,2.0668226042254458]

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
    #stats = ipopt(nlp, tol=1e-12, acceptable_tol=1e-10)
    #@time stats = ipopt(nlp)
    #@time stats = ipopt(nlp, tol=1e-2, acceptable_tol=1e-1, max_iter=1000)
    #@time stats = ipopt(nlp, tol=1e-2, max_iter=100)
    @time stats = ipopt(nlp, tol=1e-1, max_iter=1000)
    #@time stats = ipopt(nlp, tol=1e-12, max_iter=1000)
    #@time stats = ipopt(nlp, tol=1e-12, max_iter=100)
    print(stats)
    println("\nFull solution:  $(stats.solution)")

    # Note this needs to be commented out for other than the 4-node test case
    diff_solution = abs.(stats.solution - target_solution)
    println("\nTarget solution difference:  $(diff_solution)")
    println("\nMaximum solution difference:  $(findmax(diff_solution))")

    #pisolution = Pi(stats.solution)
    #println("\nPi(solution):  $(pisolution)")

    f_exp = f(target_solution)
    println("f(x_exp): $(f_exp)")

    f_sol = f(stats.solution)
    println("f(x_sol): $(f_sol)")

    for zidx in Pi_zidxs
    #  hi = h_Pi(stats.solution,zidx)
    #  println("\nh_$(zidx)(x*): $(hi)")

      hi_exp = h_Pi(target_solution,zidx)
      println("h_Pi_$(zidx)(x_exp): $(hi_exp), measurement: $(zvec[zidx])")

      #for j in keys(Ybus[zidx_nodeidx_map[zidx]])
      #  branchP = vzi(target_solution,zidx) * vj(target_solution,j)*(Gzij(zidx,j)*cos(Tzij(target_solution,zidx,j)) + Bzij(zidx,j)*sin(Tzij(target_solution,zidx,j)))
      #  i = zidx_nodeidx_map[zidx]
      #  println("\tbranch_$(i)_$(j): $(branchP)")
      #end
    end

    for zidx in Qi_zidxs
      hi_exp = h_Qi(target_solution,zidx)
      println("h_Qi_$(zidx)(x_exp): $(hi_exp), measurement: $(zvec[zidx])")
    end

    for zidx in vi_zidxs
      hi_exp = vzi(target_solution,zidx)
      println("vzi_$(zidx)(x_exp): $(hi_exp), measurement: $(zvec[zidx])")
    end

    for zidx in Ti_zidxs
      hi_exp = Tzi(target_solution,zidx)
      println("Tzi_$(zidx)(x_exp): $(hi_exp), measurement: $(zvec[zidx])")
    end
  end
#end

