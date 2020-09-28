'''
 # @ Author: Rabah Abdul Khalek
 # @ Create Time: 2020-09-28 16:33:25
 # @ Modified by: Rabah Abdul Khalek
 # @ Modified time: 2020-09-28 16:33:25
 # @ Description: Code to define max/min SFs
 '''

import os, sys, collections, yaml, lhapdf
import itertools as IT
import numpy as np
import matplotlib.pyplot as py
from pylab import *
from matplotlib import rc
from scipy.integrate import quad
cfg = lhapdf.getConfig()
cfg.set_entry("Verbosity", 0)
sys.path.append(os.getcwd()+"/../")

lhafl = {'F2NC': 908, 'FLNC': 909, 'F3NC':910}

SF_name = "NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF"
LHAPDF_dir = "/usr/local/share/LHAPDF/"

SF_path = LHAPDF_dir+SF_name

fl_to_tweak = 'FLNC'
based_on = '90CL' #options: N_SD, 68CL, 90CL.
Nrep=100

new_SF_name = SF_name+"_"+fl_to_tweak+"_"+based_on

os.system("mkdir "+new_SF_name+"_LOW")
os.system("mkdir "+new_SF_name+"_UP")

os.system("cp "+SF_path+"/"+SF_name+".info"+" "+new_SF_name+"_LOW/"+new_SF_name+"_LOW.info")
os.system("cp "+SF_path+"/"+SF_name+".info"+" "+new_SF_name+"_UP/"+new_SF_name+"_UP.info")

SFrep_list = []

for i in range(1,Nrep):

    if i<10:
        SFrep_path = SF_path+"/"+SF_name+"_000"+str(i)+".dat"
    else:
        SFrep_path = SF_path+"/"+SF_name+"_00"+str(i)+".dat"

    lines_toskip = [[1, 6], [1807, 1810], [5111, 5114], [15015, 15015]]

    acc_lines=0
    with open(SFrep_path, 'r') as f:
        lines=f
        for il,ind in enumerate(lines_toskip):
            if il==0: lines = IT.chain(IT.islice(lines, ind[0]-1), IT.islice(lines, ind[1]-(ind[0]-1), None))
            else: lines = IT.chain(IT.islice(lines, ind[0]-acc_lines-1), IT.islice(lines, ind[1]-(ind[0]-1), None))
            acc_lines += ind[1]-(ind[0]-1)

        #lines = IT.chain(IT.islice(f, 0), IT.islice(f, 6, None))
        #lines2 = IT.chain(IT.islice(lines, 1800), IT.islice(lines, 4, None))
        #lines3 = IT.chain(IT.islice(lines2, 5100), IT.islice(lines, 4, None))
        #lines4 = IT.chain(IT.islice(lines3, 15000), IT.islice(lines, 1, None))
        SFrep_data = np.genfromtxt(lines)

        SFrep_list.append(SFrep_data)

    headerf = open(SFrep_path, 'r')
    headers = headerf.readlines()
    lhaind = headers[5].split()
    lhaind = map(int, lhaind)

headerf = open(SF_path+"/"+SF_name+"_0000.dat", 'r')
headers = headerf.readlines()
lhaind = headers[5].split()
lhaind = map(int, lhaind)

if based_on == '90CL':
    SF_low = np.nanpercentile(SFrep_list, 5., axis=0)
    SF_up = np.nanpercentile(SFrep_list, 95., axis=0)
#SF_median = np.median(SFrep_list, axis=0)
SF_mean = np.mean(SFrep_list, axis=0)

#UP
os.system("cp "+SF_path+"/"+SF_name+"_0000.dat"+" "+new_SF_name+"_UP/"+new_SF_name+"_UP_0000.dat")
new_SF_file = open(new_SF_name+"_UP/"+new_SF_name+"_UP_0000.dat", "w")

max_lines=15015
max_fls = [908, 909, 910, 930, 931, 932, 940, 941, 942]
count=0

for i in range(1, max_lines+1):
    skip=False
    for ls in lines_toskip:
        ls_range = range(ls[0],ls[1]+1)
        if i in ls_range:
            skip=True
            new_SF_file.write(headers[i-1])

    if skip: continue

    for j,fl in enumerate(max_fls):
        if lhafl[fl_to_tweak] == fl:
            new_SF_file.write("{:.10e}".format(SF_up[count, j])+" ")
        else:
            new_SF_file.write(headers[i-1].split()[j]+" ")

    new_SF_file.write("\n")
    count+=1

new_SF_file.close()

#LOW
os.system("cp "+SF_path+"/"+SF_name+"_0000.dat"+" "+new_SF_name+"_LOW/"+new_SF_name+"_LOW_0000.dat")
new_SF_file = open(new_SF_name+"_LOW/"+new_SF_name+"_LOW_0000.dat", "w")

max_lines=15015
max_fls = [908, 909, 910, 930, 931, 932, 940, 941, 942]
count=0

for i in range(1, max_lines+1):
    skip=False
    for ls in lines_toskip:
        ls_range = range(ls[0],ls[1]+1)
        if i in ls_range:
            skip=True
            new_SF_file.write(headers[i-1])

    if skip: continue

    for j,fl in enumerate(max_fls):
        if lhafl[fl_to_tweak] == fl:
            new_SF_file.write("{:.10e}".format(SF_low[count, j])+" ")
        else:
            new_SF_file.write(headers[i-1].split()[j]+" ")

    new_SF_file.write("\n")
    count+=1

new_SF_file.close()

#Verify that copying was done correctly with LHAPDF
