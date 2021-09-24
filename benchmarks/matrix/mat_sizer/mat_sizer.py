#!/usr/bin/env python3

import sys
import os
import time

from csv import reader

def _main():
    with open('../mat/dimensions.csv', 'r') as read_obj:
        dims = read_obj.readline().split(',')
        xqty = int(dims[0])
        zqty = int(dims[1])

    with open('dimensions.csv', 'w') as write_obj:
        write_obj.write(str(2*xqty) + ',' + str(2*zqty)+'\n')

    with open('../mat/Supd_trip.csv', 'r') as read_obj:
        with open('Supd_trip.csv', 'w') as write_obj:
            csv_reader = reader(read_obj)
            entries = 0
            for row in csv_reader:
                i = int(row[0])
                j = int(row[1])
                val = float(row[2])
                entries += 4

                write_obj.write(str(i)+', '+str(j)+', '+str(val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j)+', '+str(2.0*val)+'\n')
                write_obj.write(str(i)+', '+str(j+zqty)+', '+str(2.0*val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j+zqty)+', '+str(3.0*val)+'\n')

            print('Supd entries: ' + str(entries))

    with open('dimensions.csv', 'a') as write_obj:
        write_obj.write(str(entries)+'\n')

    with open('../mat/K2_trip.csv', 'r') as read_obj:
        with open('K2_trip.csv', 'w') as write_obj:
            csv_reader = reader(read_obj)
            entries = 0
            for row in csv_reader:
                i = int(row[0])
                j = int(row[1])
                val = float(row[2])
                entries += 4

                write_obj.write(str(i)+', '+str(j)+', '+str(val)+'\n')
                write_obj.write(str(i+xqty)+', '+str(j)+', '+str(2.0*val)+'\n')
                write_obj.write(str(i)+', '+str(j+zqty)+', '+str(3.0*val)+'\n')
                write_obj.write(str(i+xqty)+', '+str(j+zqty)+', '+str(5.0*val)+'\n')

            print('K2 entries: ' + str(entries))

    with open('dimensions.csv', 'a') as write_obj:
        write_obj.write(str(entries)+'\n')

    with open('../mat/K3_trip.csv', 'r') as read_obj:
        with open('K3_trip.csv', 'w') as write_obj:
            csv_reader = reader(read_obj)
            entries = 0
            for row in csv_reader:
                i = int(row[0])
                j = int(row[1])
                val = float(row[2])
                entries += 4

                write_obj.write(str(i)+', '+str(j)+', '+str(val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j)+', '+str(2.0*val)+'\n')
                write_obj.write(str(i)+', '+str(j+zqty)+', '+str(3.0*val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j+zqty)+', '+str(5.0*val)+'\n')

            print('K3 entries: ' + str(entries))

    with open('dimensions.csv', 'a') as write_obj:
        write_obj.write(str(entries)+'\n')

    with open('../mat/K3sparse_trip.csv', 'r') as read_obj:
        with open('K3sparse_trip.csv', 'w') as write_obj:
            csv_reader = reader(read_obj)
            entries = 0
            for row in csv_reader:
                i = int(row[0])
                j = int(row[1])
                val = float(row[2])
                entries += 4

                write_obj.write(str(i)+', '+str(j)+', '+str(val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j)+', '+str(2.0*val)+'\n')
                write_obj.write(str(i)+', '+str(j+zqty)+', '+str(3.0*val)+'\n')
                write_obj.write(str(i+zqty)+', '+str(j+zqty)+', '+str(5.0*val)+'\n')

            print('K3sparse entries: ' + str(entries))

    with open('dimensions.csv', 'a') as write_obj:
        write_obj.write(str(entries)+'\n')


if __name__ == "__main__":
    _main()

