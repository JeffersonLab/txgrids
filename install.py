#!/usr/bin/env python
import sys, os
import argparse
from subprocess import Popen, PIPE

def get_folder_names():
    process = Popen(['ls','stf-grids'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    names=stdout.split()
    return names

def get_lhapdf_dir():
    try:
        process = Popen(['lhapdf-config','--datadir'], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stdout.strip()
    except:
        sys.exit("lhapdf-config not available")

def create_symliks():

    names=get_folder_names()
    lhapdf_dir=get_lhapdf_dir()

    cmds=['ln -s "$(pwd)/stf-grids/%s" %s/%s'%(_,lhapdf_dir,_) for _ in names]
    for cmd in cmds:
        print(cmd)
        process = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)

def remove_symliks():

    names=get_folder_names()
    lhapdf_dir=get_lhapdf_dir()

    cmds=['rm %s/%s'%(lhapdf_dir,_) for _ in names]
    for cmd in cmds:
        print(cmd)
        process = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()

if __name__=="__main__":


    ap = argparse.ArgumentParser()
    ap.add_argument('option', 
                    help='1: install symlinks 2: remove symlinks',
                    default=1,
                    type=int)
    args = ap.parse_args()

    if    args.option==1:  create_symliks()
    elif  args.option==2:  remove_symliks()
    else: sys.exit("option not available")   








