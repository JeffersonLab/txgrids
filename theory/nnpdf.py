#!/usr/bin/env python
import sys,os
import numpy as np
import params as par
import lhapdf


#--Open LHAPDF set with the structure functions (isoscalar target)
dist=lhapdf.mkPDFs("NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF")

#--Relevant EW parameters and constants
MZ = par.MZ
MW = par.MW
GF = par.GF
MN = par.MN

MZ2  = par.MZ2
MW2  = par.MW2
MN2  = par.MN2
GF2  = par.GF2
coef = GF * dist[0].quarkMass(6)**2/8/np.sqrt(2)/ np.pi**2
rho  = 1 + 3 * coef * ( 1 + coef * ( 19 - 2 * np.pi**2 ) )


#--Retrieve parameters from the LHAPDF set
Nrep = len(dist) - 1
Ql   = dist[0].q2Min**0.5
Qu   = dist[0].q2Max**0.5
xl   = dist[0].xMin
xu   = dist[0].xMax


def dsig_dx_dQ2(proc,particle,x,Q2,stot,units='pb',imem=0):
    """
    NOTE: this routine is for nu + N -> l + X,  
          need to addapt it for e + N -> l + X


    proc     = NC, CC
    particle = 12(neutrino), -22(antineutrino) 

    """

    y      = Q2 / x / ( stot - MN2 )
    omy2   = pow(1 - y, 2)
    Yplus  = 1 + omy2
    Yminus = 1 - omy2
    Q=np.sqrt(Q2)

    if proc == 'NC':

        fact=8*pow((MZ2*rho-MW2)*MW2/MZ2/Q2/rho/rho/(2-rho),2)*GF2/np.pi/x
        offset=1000

    elif proc=='CC':

        fact = GF2 / 4 / np.pi / x * pow(MW2 / ( MW2 + Q2 ), 2)
        if particle == 12:  offset=2000
        if particle ==-12:  offset=3000

    if particle== 12: nsgn= 1
    if particle==-12: nsgn=-1

    #--Conversion factor from natural unists to pb.
    if units=='pb': conv = 0.3894e9
    else: sys.exit('ERR: %s not available'%units)

    return  conv*fact*(Yplus*dist[imem].xfxQ(offset+1,x,Q)
            -y*y*dist[imem].xfxQ(offset+2,x,Q)
            +nsgn*Yminus*dist[imem].xfxQ(offset+3,x,Q))



if __name__=="__main__":


    me,Ee=0,18.0
    mN,EN=par.MN, 275.
    stot=par.get_stot(me,Ee,mN,EN)

    proc='NC'
    particle=12
    x=0.1
    imem=0
    Q2=10.0
    print(dsig_dx_dQ2(proc,particle,x,Q2,stot,units='pb',imem=0))








