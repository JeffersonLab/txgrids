#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np

import lhapdf
#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py

Q2=[5,10,100]
X1=10**np.linspace(-4,-1)
X2=np.linspace(0.101,0.99)
X=np.append(X1,X2)
iset,iF2,iFL,iF3=0,908, 909, 910 

def get_stf(X,Q2,tabname,iset,iF2,iFL,iF3):
    stf=lhapdf.mkPDF(tabname,iset)
    F2=np.array([stf.xfxQ2(iF2,x,Q2) for x in X]) 
    FL=np.array([stf.xfxQ2(iFL,x,Q2) for x in X])
    F3=np.array([stf.xfxQ2(iF3,x,Q2) for x in X])
    return F2,FL,F3


tabnames=[        
           'NNPDF31_lo_as_0118_SF'
          ,'NNPDF31_nnlo_pch_as_0118_SF'
          ,'JAM4EIC'      
         ]

data={}
for tabname in tabnames:
    data[tabname]={}
    for q2 in Q2:
        F2,FL,F3 = get_stf(X,q2,tabname,iset,iF2,iFL,iF3)
        data[tabname][q2]={}
        data[tabname][q2]['F2']=F2 
        data[tabname][q2]['FL']=FL 
        data[tabname][q2]['F3']=F3


nrows,ncols=3,3
fig = py.figure(figsize=(ncols*5,nrows*3))

cnt=0
for q2 in Q2:

    cnt+=1
    ax=py.subplot(nrows,ncols,cnt)
    for tabname in tabnames:
        F2=data[tabname][q2]['F2']
        label=tabname.replace("_","\_")
        ax.plot(X,F2,label=label)
    ax.set_ylim(0,2)
    ax.semilogx()
    if cnt==1: ax.legend(loc=2)
    if cnt==1: ax.text(0.7,0.7,r'$F_2$',size=40,transform=ax.transAxes)
    ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    if cnt==7: ax.set_xlabel(r'$x$',size=20)

    cnt+=1
    ax=py.subplot(nrows,ncols,cnt)
    for tabname in tabnames:
        FL=data[tabname][q2]['FL']
        ax.plot(X,FL)
    ax.set_ylim(0,0.3)
    ax.semilogx()
    if cnt==2: ax.text(0.7,0.7,r'$F_L$',size=40,transform=ax.transAxes)
    if cnt==8: ax.set_xlabel(r'$x$',size=20)

    cnt+=1
    ax=py.subplot(nrows,ncols,cnt)
    for tabname in tabnames:
        F3=data[tabname][q2]['F3']
        ax.plot(X,F3)#,label='F2 '+tabname,ls='-')
    ax.semilogx()
    ax.set_ylim(0,0.9)
    if cnt==3: ax.text(0.7,0.7,r'$F_3$',size=40,transform=ax.transAxes)
    if cnt==9: ax.set_xlabel(r'$x$',size=20)

    #ax.legend()
    #ax.set_ylim(0,1.7)

py.tight_layout()
py.savefig('SF_comparison.png')
py.clf()
   



