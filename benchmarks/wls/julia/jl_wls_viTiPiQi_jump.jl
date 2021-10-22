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
nlp = Model(optimizer_with_attributes(Ipopt.Optimizer,"tol"=>1e-10))

println("Start parsing files...")

nnode, nmeas, rmat, measidx_nodeidx_map, vi_measidxs, Ti_measidxs, Pi_measidxs, Qi_measidxs, YbusG, YbusB = parse_measurements()

println("Done parsing files, start defining optimization problem...")

# initialize z here then set values later when iterating over
# measurement_data.csv
zvec = Array{Float64}(undef, nmeas)

@variable(nlp,-0.5 <= v[1:nnode] <= 1.5,start=1.0)
@variable(nlp,-pi <= T[1:nnode] <= pi,start=0.0)
#@variable(nlp,v[1:nnode],start=1.0)
#@variable(nlp,T[1:nnode],start=0.0)
#@variable(nlp,-pi <= x[1:2*nnode] <= pi,start=1.0)

#@NLexpression(nlp,vi[1],sum((zvec[i] - v[measidx_nodeidx_map[i]])^2/rmat[i] for i in vi_measidxs))
#@NLexpression(nlp,Ti[1],sum((zvec[i] - T[measidx_nodeidx_map[i]])^2/rmat[i] for i in Ti_measidxs))
#@NLobjective(nlp,Min,vi[1]+Ti[1])
#print(nlp)
 
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
    zvec[imeas] = row[imeas+1]
  end
  println("*** Timestamp $(row[1]), using measurements: $(zvec)\n")

  @NLobjective(nlp, Min,
    sum((zvec[i] - v[measidx_nodeidx_map[i]])^2/rmat[i] for i in vi_measidxs) +

    sum((zvec[i] - T[measidx_nodeidx_map[i]])^2/rmat[i] for i in Ti_measidxs) +

    sum((zvec[i] - (v[measidx_nodeidx_map[i]] * sum(v[j]*(YbusG[measidx_nodeidx_map[i]][j]*cos(T[measidx_nodeidx_map[i]] - T[j]) + YbusB[measidx_nodeidx_map[i]][j]*sin(T[measidx_nodeidx_map[i]] - T[j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Pi_measidxs) +

    sum((zvec[i] - (v[measidx_nodeidx_map[i]] * sum(v[j]*(YbusG[measidx_nodeidx_map[i]][j]*sin(T[measidx_nodeidx_map[i]] - T[j]) - YbusB[measidx_nodeidx_map[i]][j]*cos(T[measidx_nodeidx_map[i]] - T[j])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Qi_measidxs))

  #@NLobjective(nlp, Min,
  #  sum((zvec[i] - x[measidx_nodeidx_map[i]])^2/rmat[i] for i in vi_measidxs) +

  #  sum((zvec[i] - x[measidx_nodeidx_map[i]+nnode])^2/rmat[i] for i in Ti_measidxs) +

  #  sum((zvec[i] - (x[measidx_nodeidx_map[i]] * sum(x[j]*(YbusG[measidx_nodeidx_map[i]][j]*cos(x[measidx_nodeidx_map[i]+nnode] - x[j+nnode]) + YbusB[measidx_nodeidx_map[i]][j]*sin(x[measidx_nodeidx_map[i]+nnode] - x[j+nnode])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Pi_measidxs) +

  #  sum((zvec[i] - (x[measidx_nodeidx_map[i]] * sum(x[j]*(YbusG[measidx_nodeidx_map[i]][j]*sin(x[measidx_nodeidx_map[i]+nnode] - x[j+nnode]) - YbusB[measidx_nodeidx_map[i]][j]*cos(x[measidx_nodeidx_map[i]+nnode] - x[j+nnode])) for j in keys(YbusG[measidx_nodeidx_map[i]]))))^2/rmat[i] for i in Qi_measidxs))

  @time optimize!(nlp)
  @time optimize!(nlp)
  solution_summary(nlp, verbose=true)
  println("v = $(value.(v))")
  println("T = $(value.(T))")
  #println("x = $(value.(x))")
end

