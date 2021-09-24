#!/usr/bin/env python3

import sys
import os
import time

import numpy as np
import scipy as sp
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import splu
from scipy.sparse.linalg import spsolve
import pymesh

from csv import reader

def process_mem_usage():
    with open('/proc/self/stat', 'r') as f:
        vsize = int(f.readline().split()[22])
        vm_usage = int(float(vsize) / (1024.0*1024.0))
        units = ' MB'
 
        return str(vm_usage) + units


def _main():
    with open('../mat/dimensions.csv', 'r') as read_obj:
        dims = read_obj.readline().split(',')
        xqty = int(dims[0])
        zqty = int(dims[1])

    Supd_lil = sp.sparse.lil_matrix((zqty, zqty))
    with open('../mat/Supd_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            Supd_lil[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Sparse Supd entries: ' + str(entries))
    Supd = Supd_lil.tocsc()

    for i in range(5):
        startTime = time.process_time()
        SupdInv = inv(Supd)
        elapsedTime = time.process_time() - startTime
        print('Python Sparse Supd inverse time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Sparse Supd inverse memory: ' + vm_used, flush=True)

    for i in range(5):
        startTime = time.process_time()
        #I = np.identity(zqty)
        #SupdInvND = splu(Supd).solve(I)
        #SupdInv = sp.sparse.csc_matrix(SupdInvND)
        #SupdInv = spsolve_triangular(SupdCSR, I)
        I = sp.sparse.eye(zqty).tocsc()
        SupdInv = spsolve(Supd, I, permc_spec='NATURAL')
        elapsedTime = time.process_time() - startTime
        print('Python Sparse Supd LU inverse time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Sparse Supd LU inverse memory: ' + vm_used, flush=True)

    for i in range(5):
        startTime = time.process_time()
        solver = pymesh.SparseSolver.create("LDLT");
        solver.compute(Supd)
        I = np.identity(zqty)
        #I = sp.sparse.eye(zqty).tocsc()
        SupdInvDense = solver.solve(I)
        #SupdInv = sp.sparse.csc_matrix(SupdInvDense)
        #print(type(SupdInv))
        elapsedTime = time.process_time() - startTime
        print('Python Sparse Supd PyMesh LDLT inverse time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Sparse Supd PyMesh LDLT inverse memory: ' + vm_used, flush=True)

    K2_lil = sp.sparse.lil_matrix((xqty, zqty))
    with open('../mat/K2_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K2_lil[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Sparse K2 entries: ' + str(entries))
    K2 = K2_lil.tocsr()
    #K2 = K2_lil.tocsc()
    #K2 = K2_lil.tocoo()
    print('Python Sparse K2 dimensions: ' + str(sp.sparse.csr_matrix.get_shape(K2)))

    startTime = time.process_time()
    K3_lil = sp.sparse.lil_matrix((zqty, zqty))
    with open('../mat/K3_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K3_lil[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Sparse K3 entries: ' + str(entries))
    elapsedTime = time.process_time() - startTime
    print('Python Sparse K3 read time: ' + str(elapsedTime), flush=True)
    K3 = K3_lil.tocsr()
    #K3 = K3_lil.tocsc()
    #K3 = K3_lil.tocoo()
    print('Python Sparse K3 dimensions: ' + str(sp.sparse.csr_matrix.get_shape(K3)))

    for i in range(5):
        startTime = time.process_time()
        #Kupd = K2 * K3
        Kupd = K2 @ K3
        elapsedTime = time.process_time() - startTime
        print('Python Sparse Kupd multiply time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Sparse Kupd multiply memory: ' + vm_used, flush=True)
    #print(Kupd)
    #print(Kupd[0,0])
    #print(Kupd[0,1])
    #print(Kupd[1,0])
    #print(Kupd[547,922])

    K3sparse_lil = sp.sparse.lil_matrix((zqty, zqty))
    with open('../mat/K3sparse_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K3sparse_lil[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Sparse K3sparse entries: ' + str(entries))
    K3sparse = K3sparse_lil.tocsr()
    print('Python Sparse K3sparse dimensions: ' + str(sp.sparse.csr_matrix.get_shape(K3sparse)))

    for i in range(5):
        startTime = time.process_time()
        Kupd_sparse = K2 @ K3sparse
        elapsedTime = time.process_time() - startTime
        print('Python Sparse Kupd_sparse multiply time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Sparse Kupd_sparse multiply memory: ' + vm_used, flush=True)


if __name__ == "__main__":
    _main()

