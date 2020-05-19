#!/usr/bin/env python
######################################
# authors: N.Sato (nsato@jlab.org) 
#
# last update: 05-19-20
#####################################
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np

import lhapdf
#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py

Q2=90.0
X1=10**np.linspace(-4,-1)
X2=np.linspace(0.101,0.99)
X=np.append(X1,X2)

def get_stf(tabname,iset,iF2,iFL,iF3):
    stf=lhapdf.mkPDF(tabname,iset)
    F2=[stf.xfxQ2(iF2,x,Q2) for x in X]  
    FL=[stf.xfxQ2(iFL,x,Q2) for x in X]
    F3=[stf.xfxQ2(iF3,x,Q2) for x in X]
    return F2,FL,F3

ax=py.subplot(111)

tabname='JAM4EIC'             
iset,iF2,iFL,iF3=0,90001,90002,90003  
F2,FL,F3 = get_stf(tabname,iset,iF2,iFL,iF3)

ax.plot(X,F2,label='F2 '+tabname,ls='-')
ax.plot(X,FL,label='FL '+tabname,ls='-')
ax.plot(X,F3,label='F3 '+tabname,ls='-')


tabname = 'NNPDF31_nnlo_pch_as_0118_SF'
iset,iF2,iFL,iF3=0,1001,1002,1003  
F2,FL,F3 = get_stf(tabname,iset,iF2,iFL,iF3)
ax.plot(X,F2,label='F2 '+tabname.replace("_","\_"),ls='--')
ax.plot(X,FL,label='FL '+tabname.replace("_","\_"),ls='--')
ax.plot(X,F3,label='F3 '+tabname.replace("_","\_"),ls='--')

ax.legend()
#ax.set_ylim(0,1.7)
ax.semilogx()
ax.set_ylabel(r'$F_x(x,Q^2=90 GeV^2)$')
ax.set_xlabel(r'$x$')

py.tight_layout()
py.savefig('SF_comparison.pdf')





