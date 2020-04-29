#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory.tools import save, load,lprint
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

wdir='.main01'


#--physical params
rs= 140.7
lum='10:fb-1'
sign=1 #--electron=1 positron=-1
ntot=10000
K=100

def veto(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1

def main00():

    #--lhapdf set and stf idx
    #tabname='JAM4EIC'             
    #iset,iF2,iFL,iF3=0,90001,90002,90003  
    tabname='NNPDF31_nnlo_pch_as_0118_rs_0.5_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   

    fname='mceg00'
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
    for k in range(K):
        lprint('%d/%d'%(k,K))
        data=mceg.gen_events(ntot)
        save(data,'%s/data00-%d.po'%(wdir,k))

def main01():

    #--lhapdf set and stf idx
    #tabname='JAM4EIC'             
    #iset,iF2,iFL,iF3=0,90001,90002,90003  
    tabname='NNPDF31_nnlo_pch_as_0118_rs_1.0_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   

    fname='mceg01'
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
    for k in range(K):
        lprint('%d/%d'%(k,K))
        data=mceg.gen_events(ntot)
        save(data,'%s/data01-%d.po'%(wdir,k))

def main_kde_test():

    ax=py.subplot(111)

    N=100000
    X=np.linspace(0,1,N)
    H,E=np.histogram(X,density=True,bins=100)
    ax.plot(E[1:],H)

    W=np.exp(-(X-0.5)**2/0.1)
    W/=np.sum(W)
    H,E=np.histogram(X,weights=W,density=True,bins=100)
    ax.plot(E[1:],H)


    kde = stats.gaussian_kde(X,weights=W)
    X=kde.resample(N)
    W=np.exp(-(X-0.5)**2/0.1)
    W/=np.sum(W)
    H,E=np.histogram(X,density=True,bins=100)
    ax.plot(E[1:],H)
    
    ax.set_ylim(0,2)
    py.savefig('test.pdf')

def main1D():
    F=os.listdir(wdir)
    data00=[load('%s/%s'%(wdir,_)) for _ in F if 'data00-' in _]
    data01=[load('%s/%s'%(wdir,_)) for _ in F if 'data01-' in _]

    d=data00[0]
    X=d['X']
    Y=d['Y']
    W=d['W']


    ax=py.subplot(111)
    H,E=np.histogram(np.log(X),weights=W,density=True,bins=100)
    ax.plot(E[1:],H)

    kde = stats.gaussian_kde(np.log(X),weights=W)
    lX=kde.resample(100000)[0]
    #I=[i for i in range(len(X)) if X[i]>0]
    H,E=np.histogram(lX,density=True,bins=100)
    ax.plot(E[1:],H)


    ax.set_xlim(-8,0)
    py.savefig('test.pdf')
    

    return 

def main2D():
    F=os.listdir(wdir)
    data00=[load('%s/%s'%(wdir,_)) for _ in F if 'data00-' in _]
    data01=[load('%s/%s'%(wdir,_)) for _ in F if 'data01-' in _]

    d=data00[0]
    X=d['X']
    Y=d['Y']
    W=d['W']

    R=((-8,0),(-7,-2))

    ax=py.subplot(111)
    #ax.hist2d(np.log(X),np.log(Y),weights=W, bins=40, norm=LogNorm(),range=R)
    H,E=np.histogram(np.log(X),weights=W,density=True,bins=100)
    ax.plot(E[1:],H,'r-')
    H,E=np.histogram(np.log(Y),weights=W,density=True,bins=100)
    ax.plot(E[1:],H,'b-')


    val = np.vstack([np.log(X),np.log(Y)])
    kde = stats.gaussian_kde(val,weights=W,bw_method=0.01)
    lX,lY=kde.resample(100000)
    H,E=np.histogram(lX,density=True,bins=100)
    ax.plot(E[1:],H,'r--')
    H,E=np.histogram(lY,density=True,bins=100)
    ax.plot(E[1:],H,'b--')

    val = np.vstack([X,Y])
    kde = stats.gaussian_kde(val,weights=W,bw_method=0.01)
    X,Y=kde.resample(100000)
    I=[i for i in range(X.size) if X[i]>0 if Y[i]>0]
    H,E=np.histogram(np.log(X[I]),density=True,bins=100)
    ax.plot(E[1:],H,'r:')
    H,E=np.histogram(np.log(Y[I]),density=True,bins=100)
    ax.plot(E[1:],H,'b:')

    ax.set_xlim(-8,0)



    #ax.hist2d(lX,lY, bins=40, norm=LogNorm(),range=R)

    #H,E=np.histogram(np.log(X),weights=W,density=True,bins=100)
    #ax.plot(E[1:],H)

    #H,E=np.histogram(np.log(X),density=True,bins=100)
    #ax.plot(E[1:],H)

    #val = np.vstack([X,Y])
    #kde = stats.gaussian_kde(val,weights=W)
    #X,Y=kde.resample(len(X)*10)
    #I=[i for i in range(len(X)) if X[i]>0  if Y[i]>0]
    #H,E=np.histogram(np.log(X[I]),density=True,bins=100)
    #ax.plot(E[1:],H)


    py.savefig('test.pdf')
    

    return 
    
    pdf=lambda x: quad(lambda y:kde([x,y]),-np.inf,np.inf)[0]
    X=np.linspace(-10,10)
    Y=[pdf(x) for x in X]
    ax=py.subplot(111)
    ax.plot(X,Y,'r-')
    ax.hist(m1,histtype='step',color='b',bins=50,normed=1)
    py.savefig('plot.pdf')


    return 

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

def main2Dhist():
    F=os.listdir(wdir)
    data00=[load('%s/%s'%(wdir,_)) for _ in F if 'data00-' in _]
    data01=[load('%s/%s'%(wdir,_)) for _ in F if 'data01-' in _]

    d=data00[0]
    X=d['X']
    Y=d['Y']
    W=d['W']

    R=((-8,0),(-7,-2))

    ax=py.subplot(121)
    ax.hist2d(np.log(X),np.log(Y),weights=W, bins=40, norm=LogNorm(),range=R)


    ax=py.subplot(122)
    val = np.vstack([np.log(X),np.log(Y)])
    kde = stats.gaussian_kde(val,weights=W,bw_method=0.01)
    lX,lY=kde.resample(10000000)
    ax.hist2d(lX,lY, bins=40, norm=LogNorm(),range=R)
    py.savefig('test.pdf')

def mainX():

    F=os.listdir(wdir)
    data00=[load('%s/%s'%(wdir,_)) for _ in F if 'data00-' in _]
    data01=[load('%s/%s'%(wdir,_)) for _ in F if 'data01-' in _]

    d0=data00[0]
    X=d0['X']
    Y=d0['Y']
    W=d0['W']

    val   = np.vstack([np.log(X),np.log(Y)])
    kde0  = stats.gaussian_kde(val,weights=W,bw_method=0.01)

    d1=data01[0]
    X=d1['X']
    Y=d1['Y']
    W=d1['W']

    val   = np.vstack([np.log(X),np.log(Y)])
    kde1  = stats.gaussian_kde(val,weights=W,bw_method=0.01)


    T0,T1=[],[]
    for _ in range(100):
        print(_)
        lXY0 = kde0.resample(1000)
        lXY1 = kde1.resample(1000)

        LL1 = np.sum(np.log([_ for _ in kde1.pdf(lXY0) if _>0]))
        LL0 = np.sum(np.log([_ for _ in kde0.pdf(lXY0) if _>0]))
        t0 = 2*(LL1-LL0)

        LL1 = np.sum(np.log([_ for _ in kde1.pdf(lXY0) if _>0]))
        LL0 = np.sum(np.log([_ for _ in kde0.pdf(lXY0) if _>0]))
        t1 = 2*(LL1-LL0)

        print(t0,t1)

        T0.append(t0)
        T1.append(t1)

    save(T0,'%s/T0A.po'%wdir)
    save(T1,'%s/T1A.po'%wdir)


def main02():


    tabname='NNPDF31_nnlo_pch_as_0118_rs_0.5_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   
    fname='mceg00'
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
    mceg0=MCEG(**data)

    tabname='NNPDF31_nnlo_pch_as_0118_rs_1.0_SF'
    iset,iF2,iFL,iF3=0,1001,1002,1003   
    fname='mceg01'
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
    mceg1=MCEG(**data)

    sig0=mceg0.get_sig_tot()
    sig1=mceg1.get_sig_tot()

    F=os.listdir(wdir)
    data00=[load('%s/%s'%(wdir,_)) for _ in F if 'data00-' in _]
    data01=[load('%s/%s'%(wdir,_)) for _ in F if 'data01-' in _]

    T0=[]
    for data in data00:
        X  = data['X']
        Y  = data['Y']
        W  = data['W']
        L0,L1=1,1
        for i in range(len(X)): 
            _L1=W[i]*mceg1._get_sigma_dxdy(X[i],Y[i],0)/sig1
            _L0=mceg0._get_sigma_dxdy(X[i],Y[i],0)/sig0
            print(_L0,_L1)
            if _L1==0: continue
            if _L0==0: continue
            if _L1>1: continue
            if _L0>1: continue
            if np.isinf(_L1): continue
            if np.isinf(_L0): continue
            if np.isnan(_L1): continue
            if np.isnan(_L0): continue
            L0*=_L0                  
            L1*=_L1
            #print(L0,L1,_L0,_L1) 
        t=2*np.log(L1/L0)
        #print t
        if np.isnan(t): continue
        T0.append(t)
    save(T0,'%s/T0.po'%wdir)

    T1=[]
    for data in data01:
        X  = data['X']
        Y  = data['Y']
        L0,L1=1,1
        for i in range(len(X)): 
            _L1=mceg1._get_sigma_dxdy(X[i],Y[i],0)/sig1
            _L0=mceg0._get_sigma_dxdy(X[i],Y[i],0)/sig0
            if _L1==0: continue
            if _L0==0: continue
            if _L1>1: continue
            if _L0>1: continue
            if np.isinf(_L1): continue
            if np.isinf(_L0): continue
            if np.isnan(_L1): continue
            if np.isnan(_L0): continue
            L0*=_L0                  
            L1*=_L1
            #print(L0,L1,_L0,_L1) 
        t=2*np.log(L1/L0)
        #print t
        if np.isnan(t): continue
        T1.append(t)
    save(T1,'%s/T1.po'%wdir)
    return 

def main03():

    T0=load('%s/T0.po'%wdir)
    T1=load('%s/T1.po'%wdir)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*5))

    ax=py.subplot(nrows,ncols,1)
    R=None
    #R=(-0.2,0)
    ax.hist(T0,bins=20,range=R,density=True,histtype='step',label='T0')
    ax.hist(T1,bins=20,range=R,density=True,histtype='step',label='T1')
    ax.legend(loc=2)
    py.tight_layout()
    py.savefig('ttest.pdf')
    return 

def main():
    ax.set_xticks(np.log([1e-4,1e-3,1e-2,1e-1]))
    ax.set_xticklabels([r'$0.0001$',r'$0.001$',r'$0.01$',r'$0.1$'])
    ax.set_yticks(np.log([1,10,100,1000,10000]))
    ax.set_yticklabels([r'$1$',r'$10$',r'$100$',r'$1000$',r'$10000$'])
    ax.set_ylabel(r'$Q^2$',size=20)
    ax.set_xlabel(r'$x$',size=20)
    ax.text(0.1,0.8,r'$\sqrt{s}=%0.2f{\rm~GeV}$'%rs,transform=ax.transAxes,size=20)

    #--make plots
    data=load('%s/data.po'%wdir)
    X  = data['X']
    Y  = data['Y']
    Q2 = data['Q2']
    W  = data['W']
    
    
    ax=py.subplot(nrows,ncols,1)
    ax.hist2d(np.log(X),np.log(Q2),weights=W, bins=40, norm=LogNorm())
    ax.set_xticks(np.log([1e-4,1e-3,1e-2,1e-1]))
    ax.set_xticklabels([r'$0.0001$',r'$0.001$',r'$0.01$',r'$0.1$'])
    ax.set_yticks(np.log([1,10,100,1000,10000]))
    ax.set_yticklabels([r'$1$',r'$10$',r'$100$',r'$1000$',r'$10000$'])
    ax.set_ylabel(r'$Q^2$',size=20)
    ax.set_xlabel(r'$x$',size=20)
    ax.text(0.1,0.8,r'$\sqrt{s}=%0.2f{\rm~GeV}$'%rs,transform=ax.transAxes,size=20)
    
    
    ax=py.subplot(nrows,ncols,2)
    ax.hist2d(np.log(X),np.log(Y),weights=W, bins=40, norm=LogNorm())
    ax.set_xticks(np.log([1e-4,1e-3,1e-2,1e-1]))
    ax.set_xticklabels([r'$0.0001$',r'$0.001$',r'$0.01$',r'$0.1$'])
    ax.set_yticks(np.log([1e-4,1e-3,1e-2,1e-1]))
    ax.set_yticklabels([r'$0.0001$',r'$0.001$',r'$0.01$',r'$0.1$'])
    ax.set_ylabel(r'$y$',size=20)
    ax.set_xlabel(r'$x$',size=20)
    
    py.tight_layout()
    py.savefig('%s/hist2d.pdf'%wdir)


if __name__=="__main__":

    #main00()
    #main01()
    #main02()
   
    #mainX()
    main03()




