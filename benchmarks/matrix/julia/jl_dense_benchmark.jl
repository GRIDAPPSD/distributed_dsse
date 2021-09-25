#!/usr/bin/env julia

function process_mem_usage()
    local vm_usage
    local units = " MB"

    open("/proc/self/stat") do io
        line = readline(io)
        vsize = parse(Float64, split(line, " ")[23]) 
        vm_usage = floor(Int, vsize / (1024.0*1024.0))
    end

    return string(vm_usage, units)
end


function dense_inv(Supd)
    SupdInv = inv(Supd)  
    return SupdInv
end


function dense_mul(K2, K3)
    Kupd = K2 * K3
    return Kupd
end


# main

open("../mat/dimensions.csv") do file
    dims = split(readline(file), ",")
    global xqty = parse(Int, dims[1])
    global zqty = parse(Int, dims[2])
end

Supd = zeros(Float64, zqty, zqty)
open("../mat/Supd_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        Supd[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end
#Supd = rand(zqty,zqty)
#println("Supd")
#println(Supd)

println("Julia Dense Supd inverse compile/gc:")
@time SupdInv = dense_inv(Supd)
@time SupdInv = dense_inv(Supd)
println("Julia Dense Supd inverse benchmark:")
@time SupdInv = dense_inv(Supd)
@time SupdInv = dense_inv(Supd)
@time SupdInv = dense_inv(Supd)
@time SupdInv = dense_inv(Supd)
@time SupdInv = dense_inv(Supd)
#println("SupdInv")
#println(SupdInv)
vm_used = process_mem_usage()
println(string("Julia Dense Supd inverse memory: ", vm_used))

K2 = zeros(Float64, xqty, zqty)
open("../mat/K2_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        K2[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end
#K2 = rand(xqty,zqty)
#println("K2")
#println(K2)

K3 = zeros(Float64, zqty, zqty)
open("../mat/K3_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        K3[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end
#K3 = rand(zqty,zqty)
#println("K3")
#println(K3)

println("Julia Dense Kupd multiply compile/gc:")
@time Kupd = dense_mul(K2,K3)
@time Kupd = dense_mul(K2,K3)
println("Julia Dense Kupd multiply benchmark:")
@time Kupd = dense_mul(K2,K3)
@time Kupd = dense_mul(K2,K3)
@time Kupd = dense_mul(K2,K3)
@time Kupd = dense_mul(K2,K3)
@time Kupd = dense_mul(K2,K3)
vm_used = process_mem_usage()
println(string("Julia Dense Kupd multiply memory: ", vm_used))
#println("Kupd")
#println(Kupd)

K3sparse = zeros(Float64, zqty, zqty)
open("../mat/K3sparse_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        K3sparse[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end

println("Julia Dense Kupd_sparse multiply benchmark:")
@time Kupd_sparse = dense_mul(K2,K3sparse)
@time Kupd_sparse = dense_mul(K2,K3sparse)
@time Kupd_sparse = dense_mul(K2,K3sparse)
@time Kupd_sparse = dense_mul(K2,K3sparse)
@time Kupd_sparse = dense_mul(K2,K3sparse)
vm_used = process_mem_usage()
println(string("Julia Dense Kupd_sparse multiply memory: ", vm_used))
#println("Kupd")

