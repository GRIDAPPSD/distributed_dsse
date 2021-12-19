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

  #nlp = ADNLPModel(f, x0)

  println("Done with defining optimization problem, start solving it...")

  # process each timestamp getting measurement data and calling solver

  # MATLAB FPI solver solution for the 4-node test case 
  #target_solution = [1.0000, 0.971576218875089, 0.956643047805180, 0.950535446461549, 0.0, 0.0, 0.0, 0.0] # FPI solver 10,000 iterations
  # MATLAB FPI solver solution for the 3p6 test case 
  #target_solution = [1.0, 1.0, 1.0, 0.9712111193225086, 0.9410836694466717, 0.9223652753951702, 0.0, -2.0943951023931957, 2.0943951023931953, -0.0375155389352858, -2.1362637169183625, 2.0491375942741321] 
  # MATLAB FPI solver solution for the 13assets_lf physical units test case 
  target_solution = [2280.1462224653337216, 2370.9304940663460002, 2291.1990472686602516, 2279.9611038436473791, 2370.9035460218810840, 2291.0555032704160112, 2362.1820781162941785, 2337.7110672203370996, 2362.1808119583511143, 2337.7095243841390584, 2277.5348581324765291, 2290.0280729393575712, 2334.7435197623253771, 2373.6298734831798356, 2340.0148421246135513, 2316.8779210656639407, 2371.6005614140126454, 2322.1821129960808321, 2269.9339329437420929, 2401.5701366959156076, 2401.6286457108431023, 2401.5788837522886752, 2270.2660272739585707, 2372.2073571348391852, 2286.7886768931412007, 2288.8924150224479490, 2330.3982535612799438, 2371.3599533499755125, 2337.0292728757217446, 2280.1464074987165986, 2370.9306694021715884, 2291.1992213259291020, 66395.3000000000029104, 66395.2999999999883585, 66395.2999999999883585, 2401.7203364564375079, 2401.7385918339664386, 2401.7263986362404466, 265.3904791881118967, 270.9028417125392707, 266.9026023653719903, -0.0480446300322439, -2.1104985943393637, 2.0545258974007106, -0.0480650512861194, -2.1104993103261331, 2.0545086644810366, -2.1086030027029521, 2.0722242874676491, -2.1086034364935085, 2.0722245481961461, -0.0480391337561771, 2.0529656575179640, -0.0229219642187778, -2.1069033570263000, 2.0719690478578934, -0.0311579982700007, -2.1080452706239958, 2.0654498954911631, -0.0476869700719445, -0.0000394796490696, -2.0944289185520728, 2.0943476493206190, -0.0494099448676698, -2.1112820332230458, 2.0554359284967827, 2.0509730331871685, -0.0234745221225297, -2.1072314747030054, 2.0717780313744405, -0.0480446535613301, -2.1104986233690131, 2.0545258706446399, 0.5235987755982988, -1.5707963267948963, 2.6179938779914944, -0.0000273674186166, -2.0944166624406009, 2.0943655287652119, -0.0297801199764095, -2.1114943921009415, 2.0673840447661922]

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
    #@time stats = ipopt(nlp, tol=1e-10)
    #@time stats = ipopt(nlp, tol=1e-12, max_iter=1000)
    #@time stats = ipopt(nlp, tol=1e-12, max_iter=100)
    #print(stats)
    #println("\nFull solution:  $(stats.solution)")

    # Note this needs to be commented out for other than the 4-node test case
    #diff_solution = stats.solution - target_solution
    #println("\nTarget solution difference:  $(diff_solution)")

    #pisolution = Pi(stats.solution)
    #println("\nPi(solution):  $(pisolution)")

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

