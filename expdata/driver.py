#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py
from theory.tools import checkdir,save,load
import pandas as pd
from scipy.interpolate import griddata

def gen_grid(path=""):
    L=open(path+'src/xQ2binTable-xiaoxuan-060220.dat').readlines()
    L=[_.strip() for _ in L]
    L=[_ for _ in L if _!='']
    #for _ in L: print(_)    

    X   =np.array([float(_) for _ in L[1].split()])
    Q2  =np.array([float(_) for _ in L[3].split()])
    X0  =np.exp(0.5*(np.log(X[1:]) +np.log(X[:-1])))
    Q20 =np.exp(0.5*(np.log(Q2[1:])+np.log(Q2[:-1])))

    #--gen grid
    rs    = 140.7
    s     = rs**2
    W2min = 10.0    
    Q2min = 1
    xmin  = Q2min/s
    M     = 0.938

    grid=[]
    for x in X:
        for q2 in Q2:
            if q2<Q2min: continue
            if x<xmin: continue
            Q2max=s*x
            if q2>Q2max: continue
            W2=M**2 + q2*(1-x)/x
            if W2<W2min: continue
            #--bin passes the filters
            grid.append([x,q2])

    grid0=[]
    for x in X0:
        for q2 in Q20:
            if q2<Q2min: continue
            if x<xmin: continue
            Q2max=s*x
            if q2>Q2max: continue
            W2=M**2 + q2*(1-x)/x
            if W2<W2min: continue
            #--bin passes the filters
            grid0.append([x,q2])

    return grid0, grid

def main00():

    grid0, grid=gen_grid()
 
    BINS=[] 
    for bin0 in grid0:
        x0,q20  = bin0
        coords  = []
        R2      = []
        for kin in grid:
            x,q2=kin
            r2=(np.log(x)-np.log(x0))**2+(np.log(q2)-np.log(q20))**2
            coords.append(kin)
            R2.append(r2)
        I=np.argsort(R2)
        if len(I)<4: continue
        #print([R2[i] for i in I[:4]])
        rel1=np.abs((R2[I[0]]-R2[I[1]])/R2[I[1]])*100
        rel2=np.abs((R2[I[2]]-R2[I[3]])/R2[I[2]])*100
        if rel1>10: continue
        if rel2>10: continue
        coords=[coords[i] for i in I[:4]]
        X,Q2 =np.transpose(coords)
        if len(np.unique(X))!=2: continue
        if len(np.unique(Q2))!=2: continue
    
        BINS.append(coords)

    checkdir('data')
    np.save('data/xQ2binTable-xiaoxuan-060220.npy',BINS)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*4))
    ax=py.subplot(nrows,ncols,1)

    for kin in grid:
        x,q2=kin
        ax.plot(x,q2,'k.')

    for kin in grid0:
        x,q2=kin
        ax.plot(x,q2,'r.')

    for bin in BINS:
        x,q2 =np.transpose(bin)
        xmin =np.amin(x)
        xmax =np.amax(x)
        q2min=np.amin(q2)
        q2max=np.amax(q2)
        ax.plot([xmin,xmax,xmax,xmin,xmin],[q2min,q2min,q2max,q2max,q2min])


    ax.set_ylabel(r'$Q^2$',size=20)
    ax.set_xlabel(r'$x_{\rm bj}$',size=20)
    ax.semilogx()
    ax.semilogy()
    py.tight_layout()
    checkdir('gallery')
    py.savefig('gallery/main00.pdf')

def main01():

    scale=300

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*4))
    ax=py.subplot(nrows,ncols,1)

    fname = 'tables/expdata/10026.xlsx'
    tab   = pd.read_excel(fname)
    tab   = tab.to_dict(orient='list')
    val   = np.array(tab['value'])
    stat  = np.array(tab['stat_u'])
    syst  = np.array(tab['syst_u'])
    alpha = np.sqrt(syst**2)
    LX0   = np.log(tab['X'])
    LQ20  = np.log(tab['Q2'])
    rel=alpha/val
    for i in range(len(tab['X'])):
        ax.plot(tab['X'][i],tab['Q2'][i],'ko',alpha=0.2,markersize=scale*rel[i],mfc='None')


    grid0, grid=gen_grid()
    rel_list=[]
    data=[]
    for kin in grid0:
        x,q2=kin
        alpha_bin=griddata((LX0,LQ20),alpha,(np.log(x) ,np.log(q2)), fill_value=0, method='cubic')
        val_bin=griddata((LX0,LQ20),val,(np.log(x) ,np.log(q2)), fill_value=0, method='cubic')
        rel=alpha_bin/val_bin
        if np.isnan(rel)==False:
            ax.plot(x,q2,'ro',markersize=scale*rel,mfc='None')
            rel_list.append(rel)
            data.append([x,q2,rel])

    for kin in grid0:
        x,q2=kin
        alpha_bin=griddata((LX0,LQ20),alpha,(np.log(x) ,np.log(q2)), fill_value=0, method='cubic')
        val_bin=griddata((LX0,LQ20),val,(np.log(x) ,np.log(q2)), fill_value=0, method='cubic')
        rel=alpha_bin/val_bin
        if np.isnan(rel)==True:
            rel=np.amax(rel_list)
            ax.plot(x,q2,'bo',markersize=scale*rel,mfc='None')
            data.append([x,q2,rel])

    checkdir('data')
    np.save('data/xQ2binTable-xiaoxuan-060220+syst.npy',data)

    ax.set_ylabel(r'$Q^2$',size=20)
    ax.set_xlabel(r'$x_{\rm bj}$',size=20)
    ax.semilogx()
    ax.semilogy()
    py.tight_layout()
    checkdir('gallery')
    py.savefig('gallery/main01.pdf')


if __name__== "__main__":

    #main00()
    main01()




##-deprecated
#fname = 'tables/expdata/10016.xlsx'
#tab   = pd.read_excel(fname)
#tab   = tab.to_dict(orient='list')
#val   = np.array(tab['value'])
#stat  = np.array(tab['stat_u'])
#syst  = np.array(tab['syst_u'])
#alpha = np.sqrt(stat**2+syst**2)
#LX0   = np.log(tab['X'])
#LQ20  = np.log(tab['Q2'])
#ax.plot(tab['X'],tab['Q2'],'k.')

#fname = 'tables/expdata/10010.xlsx'
#tab   = pd.read_excel(fname)
#tab   = tab.to_dict(orient='list')
#val   = np.array(tab['value'])
#stat  = np.array(tab['stat_u'])
#syst  = np.array(tab['syst_u'])
#alpha = np.sqrt(stat**2+syst**2)
#LX0   = np.log(tab['X'])
#LQ20  = np.log(tab['Q2'])
#ax.plot(tab['X'],tab['Q2'],'k.')


#fname = 'tables/expdata/10020.xlsx'
#tab   = pd.read_excel(fname)
#tab   = tab.to_dict(orient='list')
#val   = np.array(tab['value'])
#stat  = np.array(tab['stat_u'])
#syst  = np.array(tab['syst_u'])
#alpha = np.sqrt(stat**2+syst**2)
#LX0   = np.log(tab['X'])
#LQ20  = np.log(tab['Q2'])
#ax.plot(tab['X'],tab['Q2'],'k.')
