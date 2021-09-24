#!/usr/bin/env python3

import sys
import os
import time

import scipy as sp
from scipy.sparse.linalg import inv

from csv import reader

def _main():
    K3_rand = sp.sparse.random(923, 923, density=0.02)

    with open('K3sparse_trip.csv', 'w') as write_obj:
        for i,j,v in zip(K3_rand.row, K3_rand.col, K3_rand.data):
            write_obj.write(str(i)+', '+str(j)+', '+str(v)+'\n')


if __name__ == "__main__":
    _main()

