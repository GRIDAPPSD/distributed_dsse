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

SupdSparse = sparse(Supd)

SupdCH = SuiteSparse.CHOLMOD.Sparse(SupdSparse);

#=
println("Julia Sparse Supd inverse compile/gc:")
@time SupdInv = sparse_inv(SupdCH, zqty)
@time SupdInv = sparse_inv(SupdCH, zqty)
println("Julia Sparse Supd inverse benchmark:")
@time SupdInv = sparse_inv(SupdCH, zqty)
@time SupdInv = sparse_inv(SupdCH, zqty)
@time SupdInv = sparse_inv(SupdCH, zqty)
@time SupdInv = sparse_inv(SupdCH, zqty)
@time SupdInv = sparse_inv(SupdCH, zqty)
vm_used = process_mem_usage()
println(string("Julia SuiteSparse Supd inverse memory: ", vm_used))
#println("SupdInv")
#println(SupdInv)
=#
println("*******************************************************")
println("*** SuiteSparse Supd inverse not supported in Julia ***")
println("*******************************************************")

println("Julia Sparse Supd dense inverse compile/gc:")
@time SupdInv = sparse_inv_dense(SupdCH)
@time SupdInv = sparse_inv_dense(SupdCH)
println("Julia Sparse Supd dense inverse benchmark:")
@time SupdInv = sparse_inv_dense(SupdCH)
@time SupdInv = sparse_inv_dense(SupdCH)
@time SupdInv = sparse_inv_dense(SupdCH)
@time SupdInv = sparse_inv_dense(SupdCH)
@time SupdInv = sparse_inv_dense(SupdCH)
vm_used = process_mem_usage()
println(string("Julia SuiteSparse Supd inverse dense memory: ", vm_used))

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

K2Sparse = sparse(K2)
K3Sparse = sparse(K3)

K2CH = SuiteSparse.CHOLMOD.Sparse(K2Sparse);
K3CH = SuiteSparse.CHOLMOD.Sparse(K3Sparse);

println("Julia SuiteSparse Kupd multiply compile/gc:")
@time KupdCH = sparse_mul(K2CH,K3CH)
@time KupdCH = sparse_mul(K2CH,K3CH)
println("Julia SuiteSparse Kupd multiply benchmark:")
@time KupdCH = sparse_mul(K2CH,K3CH)
@time KupdCH = sparse_mul(K2CH,K3CH)
@time KupdCH = sparse_mul(K2CH,K3CH)
@time KupdCH = sparse_mul(K2CH,K3CH)
@time KupdCH = sparse_mul(K2CH,K3CH)
vm_used = process_mem_usage()
println(string("Julia SuiteSparse Kupd multiply memory: ", vm_used))
#println("KupdCH")
#println(KupdCH)

K3sparse = zeros(Float64, zqty, zqty)
open("../mat/K3sparse_trip.csv") do file
    for li in eachline(file)
        trip = split.(li, ",")
        K3sparse[parse(Int,trip[1])+1,parse(Int,trip[2])+1] = parse(Float64,trip[3])
    end
end

K3ReallySparse = sparse(K3sparse)
K3ReallyCH = SuiteSparse.CHOLMOD.Sparse(K3ReallySparse);

println("Julia SuiteSparse Kupd_sparse multiply benchmark:")
@time KupdCH_sparse = sparse_mul(K2CH,K3ReallyCH)
@time KupdCH_sparse = sparse_mul(K2CH,K3ReallyCH)
@time KupdCH_sparse = sparse_mul(K2CH,K3ReallyCH)
@time KupdCH_sparse = sparse_mul(K2CH,K3ReallyCH)
@time KupdCH_sparse = sparse_mul(K2CH,K3ReallyCH)
vm_used = process_mem_usage()
println(string("Julia SuiteSparse Kupd_sparse multiply memory: ", vm_used))

