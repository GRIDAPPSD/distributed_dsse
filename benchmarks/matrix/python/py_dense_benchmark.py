#!/usr/bin/env python3

import sys
import os
import time

import numpy as np

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

    Supd = np.zeros((zqty, zqty))
    with open('../mat/Supd_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            Supd[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Dense Supd entries: ' + str(entries))

    for i in range(5):
        startTime = time.process_time()
        SupdInv = np.linalg.inv(Supd)
        elapsedTime = time.process_time() - startTime
        print('Python Dense Supd inverse time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Dense Supd inverse memory: ' + vm_used, flush=True)

    K2 = np.zeros((xqty, zqty))
    with open('../mat/K2_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K2[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Dense K2 entries: ' + str(entries))

    startTime = time.process_time()
    K3 = np.zeros((zqty, zqty))
    with open('../mat/K3_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K3[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Dense K3 entries: ' + str(entries))
    elapsedTime = time.process_time() - startTime
    print('Python Dense K3 read time: ' + str(elapsedTime), flush=True)

    for i in range(5):
        startTime = time.process_time()
        Kupd = np.matmul(K2, K3)
        elapsedTime = time.process_time() - startTime
        print('Python Dense Kupd multiply time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Dense Kupd multiply memory: ' + vm_used, flush=True)
    #print(Kupd)
    #print(Kupd[0,0])
    #print(Kupd[0,1])
    #print(Kupd[1,0])
    #print(Kupd[547,922])

    K3_sparse = np.zeros((zqty, zqty))
    with open('../mat/K3sparse_trip.csv', 'r') as read_obj:
        csv_reader = reader(read_obj)
        entries = 0
        for row in csv_reader:
            K3_sparse[int(row[0]),int(row[1])] = float(row[2])
            entries += 1
        print('Python Dense K3sparse entries: ' + str(entries))

    for i in range(5):
        startTime = time.process_time()
        Kupd_sparse = np.matmul(K2, K3_sparse)
        elapsedTime = time.process_time() - startTime
        print('Python Dense Kupd_sparse multiply time: ' + str(elapsedTime), flush=True)
    vm_used = process_mem_usage()
    print('Python Dense Kupd_sparse multiply memory: ' + vm_used, flush=True)


if __name__ == "__main__":
    _main()

