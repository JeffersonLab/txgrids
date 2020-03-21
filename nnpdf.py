#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
import lhapdf

#--Open LHAPDF set with the structure functions (isoscalar target)
dist=lhapdf.mkPDFs("NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF")

#--Relevant EW parameters and constants
MZ   = 91.1876
MW   = 80.385
GF   = 1.1663787e-5
MN   = 0.938272046

#--Derived parameters.See Eqs. (8.22) and (8.25) of
#--https://arxiv.org/pdf/0709.1075.pdfhttps://arxiv.org/pdf/0709.1075.pdf
MZ2  = MZ * MZ
MW2  = MW * MW
MN2  = MN * MN
GF2  = GF * GF
coef = GF * dist[0].quarkMass(6)**2/8/np.sqrt(2)/ np.pi**2
rho  = 1 + 3 * coef * ( 1 + coef * ( 19 - 2 * np.pi**2 ) )


#--Retrieve parameters from the LHAPDF set
Nrep = len(dist) - 1
Ql   = dist[0].q2Min**0.5
Qu   = dist[0].q2Max**0.5
xl   = dist[0].xMin
xu   = dist[0].xMax

Enu=30.0
stot = MN2 + 2 * MN * Enu

#--Conversion factor from natural unists to pb.
conv = 0.3894e9


def dsigdlnxdlnQ2(proc,particle,x,Q2,stot,imem):
    """
    proc     = NC, CC
    particle = 12(neutrino), -22(antineutrino) 
    """

    y      = Q2 / x / ( stot - MN2 )
    omy2   = pow(1 - y, 2)
    Yplus  = 1 + omy2
    Yminus = 1 - omy2
    Q=np.sqrt(Q2)

    if proc == 'NC':

        fact=conv*8*pow((MZ2*rho-MW2)*MW2/MZ2/Q2/rho/rho/(2-rho),2)*GF2/np.pi/x
        offset=1000

    elif proc=='CC':

        fact = conv * GF2 / 4 / np.pi / x * pow(MW2 / ( MW2 + Q2 ), 2)
        if particle == 12:  offset=2000
        if particle ==-12:  offset=3000

    if particle== 12: nsgn= 1
    if particle==-12: nsgn=-1


    #--The factor x is the jacobian of the x
    #--integration due to the fact that dx = x *
    #--d(log(x)).
    return  x*fact*(Yplus*dist[imem].xfxQ(offset+1,x,Q)
            -y*y*dist[imem].xfxQ(offset+2,x,Q)
            +nsgn*Yminus*dist[imem].xfxQ(offset+3,x,Q))


proc='NC'
particle=12
x=0.1
imem=0
Q2=10.0


print(dsigdlnxdlnQ2(proc,particle,x,Q2,stot,imem))







#if __name__=="__main__":
#
#    name="NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF"
#    #st=load_stfuncs(name)









