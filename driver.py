#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import theory as thy
from theory import params as par
from theory import nnpdf

def main():

    #--setup beams energies
    me,Ee=0,18.0
    mN,EN=par.MN, 275.

    #--get s=(pA+pB)^2
    stot=par.get_stot(me,Ee,mN,EN)

    proc='NC'
    particle=12  #--neutrino

    x=0.1
    imem=0
    Q2=10.0

    print(thy.nnpdf.dsig_dx_dQ2(proc,particle,x,Q2,stot,units='pb',imem=0))


if __name__=="__main__":

    main()




