#!/usr/bin/env julia

using SparseArrays, SuiteSparse, LinearAlgebra

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


function sparse_inv(Supd, zqty)
    # Julia recognizes this as trying to compute the inverse and errors out
    #SupdInv = Supd \ I 

    eye = Matrix{Float64}(I, zqty, zqty)
    #SupdInvDense = Supd \ eye
    SupdSym = Symmetric(Supd)
    SupdInvDense = SupdSym \ eye
    SupdInv = sparse(SupdInvDense)
    return SupdInv
end


function sparse_inv_dense(Supd)
    SupdDense = Matrix(Supd)
    SupdInvDense = inv(SupdDense)
    SupdInv = sparse(SupdInvDense)
    return SupdInv
end


function sparse_mul(K2, K3)
    Kupd = K2 * K3
    return Kupd
end


# main

open("../mat/dimensions.csv") do file
    dims = split(readline(file), ",")
    global xqty = parse(Int, dims[1])
    global zqty = parse(Int, dims[2])
    global Supd_entries = parse(Int, readline(file))
    global K2_entries = parse(Int, readline(file))
    global K3_entries = parse(Int, readline(file))
    global K3sparse_entries = parse(Int, readline(file))
end


Supd = spzeros(Float64, zqty, zqty)
println("Julia Sparse reading Supd...")
open("../mat/Supd_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        Supd[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end
#Supd = rand(zqty,zqty)
#println("Supd")
#println(Supd)

println("Julia Sparse Supd inverse compile/gc:")
@time SupdInv = sparse_inv(Supd, zqty)
@time SupdInv = sparse_inv(Supd, zqty)
println("Julia Sparse Supd inverse benchmark:")
@time SupdInv = sparse_inv(Supd, zqty)
@time SupdInv = sparse_inv(Supd, zqty)
@time SupdInv = sparse_inv(Supd, zqty)
@time SupdInv = sparse_inv(Supd, zqty)
@time SupdInv = sparse_inv(Supd, zqty)
vm_used = process_mem_usage()
println(string("Julia Sparse Supd inverse memory: ", vm_used))
#println("SupdInv")
#println(SupdInv)

println("Julia Sparse Supd dense inverse compile/gc:")
@time SupdInv = sparse_inv_dense(Supd)
@time SupdInv = sparse_inv_dense(Supd)
println("Julia Sparse Supd dense inverse benchmark:")
@time SupdInv = sparse_inv_dense(Supd)
@time SupdInv = sparse_inv_dense(Supd)
@time SupdInv = sparse_inv_dense(Supd)
@time SupdInv = sparse_inv_dense(Supd)
@time SupdInv = sparse_inv_dense(Supd)
vm_used = process_mem_usage()
println(string("Julia Sparse Supd dense inverse memory: ", vm_used))

K2 = zeros(Float64, xqty, zqty)
println("Julia Sparse reading K2...")
open("../mat/K2_trip.csv") do file
    local entries = 0
    for li in eachline(file)
        trip = split.(li, ",")
        K2[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
        entries += 1
        if mod(entries, 10000) == 0
            perdone = convert(Int, floor(100.0*(entries/K2_entries)))
            println(string("Julia Sparse read ",entries," K2 entries (",perdone,"%)"))
        end
    end
end
#K2 = rand(xqty,zqty)
#println("K2")
#println(K2)

println("Julia Sparse reading K3...")
K3 = zeros(Float64, zqty, zqty)
open("../mat/K3_trip.csv") do file
    local entries = 0
    for li in eachline(file)
        trip = split.(li, ",")
        K3[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
        entries += 1
        if mod(entries, 10000) == 0
            perdone = convert(Int, floor(100.0*(entries/K3_entries)))
            println(string("Julia Sparse read ",entries," K3 entries (",perdone,"%)"))
        end
    end
end
#K3 = rand(zqty,zqty)
#println("K3")
#println(K3)

println("Julia Sparse converting K2 to sparse matrix...")
K2Sparse = sparse(K2)
println("Julia Sparse converting K3 to sparse matrix...")
K3Sparse = sparse(K3)

println("Julia Sparse Kupd multiply compile/gc:")
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
println("Julia Sparse Kupd multiply benchmark:")
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
@time Kupd = sparse_mul(K2Sparse,K3Sparse)
vm_used = process_mem_usage()
println(string("Julia Sparse Kupd multiply memory: ", vm_used))
#println("Kupd")
#println(Kupd)

println("Julia Sparse reading K3sparse...")
K3sparse = zeros(Float64, zqty, zqty)
open("../mat/K3sparse_trip.csv") do file
    local entries = 0
    for li in eachline(file)
        trip = split.(li, ",")
        K3sparse[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
        entries += 1
        if mod(entries, 10000) == 0
            perdone = convert(Int, floor(100.0*(entries/K3sparse_entries)))
            println(string("Julia Sparse read ",entries," K3sparse entries (",perdone,"%)"))
        end
    end
end

println("Julia Sparse converting K3sparse to sparse matrix...")
K3ReallySparse = sparse(K3sparse)

println("Julia Sparse Kupd_sparse multiply benchmark:")
@time Kupd_sparse = sparse_mul(K2Sparse,K3ReallySparse)
@time Kupd_sparse = sparse_mul(K2Sparse,K3ReallySparse)
@time Kupd_sparse = sparse_mul(K2Sparse,K3ReallySparse)
@time Kupd_sparse = sparse_mul(K2Sparse,K3ReallySparse)
@time Kupd_sparse = sparse_mul(K2Sparse,K3ReallySparse)
vm_used = process_mem_usage()
println(string("Julia Sparse Kupd_sparse multiply memory: ", vm_used))

