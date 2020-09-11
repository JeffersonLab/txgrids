#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory.tools import save, load,lprint,checkdir
from theory.mceg  import MCEG

#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

from scipy import stats
from scipy.integrate import quad
import gkde

wdir='.ttest'

#--physical params
rs= 140.7
lum='10:fb-1'
sign=1 #--electron=1 positron=-1
ntot=10000

def veto(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1

def gen_events():

    #--null hypotesis
    tabname='NNPDF31_nnlo_pch_as_0118_rs_1.0_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   
    fname='mceg0'
    data={}
    data['wdir']    =  wdir   
    data['tabname'] =  tabname
    data['iset']    =  iset   
    data['iF2']     =  iF2    
    data['iFL']     =  iFL    
    data['iF3']     =  iF3    
    data['sign']    =  sign   
    data['rs']      =  rs     
    data['fname']   =  fname  
    data['veto']    =  veto
    mceg=MCEG(**data)
    mceg.buil_mceg()
    data=mceg.gen_events(ntot)
    save(data,'%s/evt0.po'%(wdir))

    #--alternative hypotesis
    tabname='NNPDF31_nnlo_pch_as_0118_rs_0.5_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   
    fname='mceg1'
    data={}
    data['wdir']    =  wdir   
    data['tabname'] =  tabname
    data['iset']    =  iset   
    data['iF2']     =  iF2    
    data['iFL']     =  iFL    
    data['iF3']     =  iF3    
    data['sign']    =  sign   
    data['rs']      =  rs     
    data['fname']   =  fname  
    data['veto']    =  veto
    mceg=MCEG(**data)
    mceg.buil_mceg()
    data=mceg.gen_events(ntot)
    save(data,'%s/evt1.po'%(wdir))

def gen_tvals():

    #--retrieve events
    evt0=load('%s/evt0.po'%wdir)
    evt1=load('%s/evt1.po'%wdir)

    #--build kde for evt0
    X=evt0['X']
    Y=evt0['Y']
    W=evt0['W']
    val   = np.vstack([np.log(X),np.log(Y)])
    kde0  = stats.gaussian_kde(val,weights=W,bw_method=0.01)

    #--build kde for evt1
    X=evt1['X']
    Y=evt1['Y']
    W=evt1['W']
    val   = np.vstack([np.log(X),np.log(Y)])
    kde1  = stats.gaussian_kde(val,weights=W,bw_method=0.01)

    #--gen t0 and t1
    T0,T1 = [],[]
    for _ in range(100):
        print(_)
        lXY0 = kde0.resample(1000)
        lXY1 = kde1.resample(1000)

        LL1 = np.sum(kde1.logpdf(lXY0))
        LL0 = np.sum(kde0.logpdf(lXY0))
        t0 = 2*(LL1-LL0)

        LL1 = np.sum(kde1.logpdf(lXY1))
        LL0 = np.sum(kde0.logpdf(lXY1))
        t1 = 2*(LL1-LL0)

        print(t0,t1)

        T0.append(t0)
        T1.append(t1)

    save(T0,'%s/T0.po'%wdir)
    save(T1,'%s/T1.po'%wdir)

def plot_tvals():

    T0=load('%s/T0.po'%wdir)
    T1=load('%s/T1.po'%wdir)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*5))

    ax=py.subplot(nrows,ncols,1)
    R=None
    #R=(-0.2,0)
    ax.hist(T0,bins=10,range=R,density=True,histtype='step',label='T0')
    ax.hist(T1,bins=10,range=R,density=True,histtype='step',label='T1')
    ax.set_xlabel(r'$t$',size=20)
    ax.set_ylabel(r'$\rm Normalized~Yield$',size=20)
    ax.legend(loc=1)
    py.tight_layout()
    py.savefig('%s/ttest.pdf'%(wdir))


if __name__=="__main__":

    #gen_events()
    #gen_tvals()
    plot_tvals()



