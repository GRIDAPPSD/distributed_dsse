#!/usr/bin/env python3

import math

def noground(x):
  return x

def ground(x):
  return round(x*100.0)/100.0

# read vnom.csv into dictionary
#Nodename,Mag,Arg
#SOURCEBUS.1,66395.3,30.0
#SOURCEBUS.2,66395.3,-90.0
#SOURCEBUS.3,66395.3,150.0
#671.1,2401.65,0.0
#671.2,2401.63,-120.0
#671.3,2401.61,120.0

VnomMag = {}
VnomArg = {}

with open('vnom.csv') as vnomfile:
  for line in vnomfile:
    fields = line.split(',')
    if fields[0] == 'Nodename':
      continue
    VnomMag[fields[0]] = float(fields[1])
    VnomArg[fields[0]] = math.radians(float(fields[2]))
#print(VnomMag)
#print(VnomArg)

# read and convert measurements.csv file
#ztype,zid,znode1,znode2,zval-change,zsig-change,zpseudo,znomval-change
#Pi,pseudo_P_671.1,671.1,671.1,-0.1925,0.194613,1,-0.1925
#Qi,pseudo_Q_671.1,671.1,671.1,-0.11,0.111282,1,-0.11

sbase = 1.0e+6
idx = 1
MeasType = {}
MeasNode = {}

newmeasfile = open('measurements.csv.new', 'w')

with open('measurements.csv') as measfile:
  for line in measfile:
    fields = line.split(',')
    if fields[0] == 'ztype':
      newmeasfile.write(line)
      continue

    MeasType[idx] = fields[0]
    MeasNode[idx] = fields[2]

    if fields[0] == 'vi':
      zval = ground(float(fields[4]) * VnomMag[fields[2]])
      zsig = ground(float(fields[5]) * VnomMag[fields[2]])
      znomval = ground(float(fields[7]) * VnomMag[fields[2]])
    elif fields[0] == 'Ti':
      zval = ground(float(fields[4]))
      zsig = ground(float(fields[5]))
      znomval = ground(float(fields[7]))
      #zval = ground(float(fields[4]) + VnomArg[fields[2]])
      #zsig = ground(float(fields[5]) + VnomArg[fields[2]])
      #znomval = ground(float(fields[7]) + VnomArg[fields[2]])
    else: # Pi and Qi
      zval = ground(float(fields[4]) * sbase)
      zsig = ground(float(fields[5]) * sbase)
      znomval = ground(float(fields[7]) * sbase)

    newmeasfile.write(fields[0]+','+fields[1]+','+fields[2]+','+fields[3]+','+str(zval)+','+str(zsig)+','+fields[6]+','+str(znomval)+'\n')

    idx += 1

newmeasfile.close()

newdatafile = open('measurement_data.csv.new', 'w')

with open('measurement_data.csv') as datafile:
  for line in datafile:
    fields = line.split(',')
    if fields[0] == 'timestamp':
      newdatafile.write(line)
      continue

    newdatafile.write(fields[0])

    for idx in range(1, len(fields)):
      if MeasType[idx] == 'vi':
        zval = ground(float(fields[idx]) * VnomMag[MeasNode[idx]])
      elif MeasType[idx] == 'Ti':
        zval = ground(float(fields[idx]) + VnomArg[MeasNode[idx]])
      else: # Pi or Qi
        zval = ground(float(fields[idx]) * sbase)

      newdatafile.write(','+str(zval))
    newdatafile.write('\n')

newdatafile.close()


