#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
import pandas as pd
#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

from theory.tools import checkdir,save,load
from theory.idis  import IDIS

tabname='JAM4EIC'             

def veto(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1

def get_hera_bins():
    xbins=[]
    xbins.append([3.5e-5  , 5.5e-5   ])
    xbins.append([5.5e-5  , 9.5e-5   ])
    xbins.append([9.5e-5  , 1.5e-4   ])
    xbins.append([1.5e-4  , 2.5e-4   ])
    xbins.append([2.5e-4  , 3.5e-4   ])
    xbins.append([3.5e-4  , 5.5e-4   ])
    xbins.append([5.5e-4  , 6.6e-4   ])
    xbins.append([7.5e-4  , 8.5e-4   ])
    xbins.append([8.5e-4  , 0.98e-3  ])
    xbins.append([0.98e-3 , 1.9e-3   ])
    xbins.append([0.0019  , 0.0030   ])
    xbins.append([0.0030  , 0.0040   ])
    xbins.append([4.8e-3  , 5.6e-3   ])
    xbins.append([7.9e-3  , 8.6e-3   ])
    xbins.append([1.0e-2  , 1.5e-2   ])
    xbins.append([1.9e-2  , 2.2e-2   ])
    xbins.append([3.0e-2  , 3.4e-2   ])
    xbins.append([4.8e-2  , 5.6e-2   ])
    xbins.append([7.9e-2  , 8.8e-2   ])
    xbins.append([12e-2   , 14e-2    ])
    xbins.append([17e-2   , 19e-2    ])
    xbins.append([23e-2   , 26e-2    ])
    xbins.append([38.0e-2 , 42.0e-2  ])
    xbins.append([62.0e-2 , 66e-2    ])
    return xbins[::-1]

data={}
data['tabname'] = tabname
data['iset']    = 0 
data['iF2']     = 908
data['iFL']     = 909    
data['iF3']     = 910    
data['sign']    = 1 #--electron=1 positron=-1
data['veto']    = veto

idis=IDIS(**data)
rs=320.0
exp=pd.read_excel('expdata/hera.xlsx')
xbins=get_hera_bins()    

nrows,ncols=1,1
fig = py.figure(figsize=(ncols*4,nrows*5))

ax=py.subplot(nrows,ncols,1)
for ibin in range(len(xbins)):

    factor=2**ibin

    xmin,xmax=xbins[ibin]
    xbin=exp.query('%f<X and X<%f and Q2>1'%(xmin,xmax))
    alpha=np.sqrt(xbin['stat_u']**2+xbin['syst_u']**2)
    hhera=ax.errorbar(xbin.Q2,xbin.value*factor
                     ,alpha*factor,fmt='k.',markersize=2,elinewidth=1)
    
    thy=np.array([idis.get_sigma_red(xbin.X.values[i],xbin.Q2.values[i],rs,-1) 
                  for i in range(len(xbin.X.values))])

    hthy,=ax.plot(xbin.Q2,thy*factor,ls='-')
  

#ax.set_ylim(1e2,1e7)
ax.set_xlim(1, 5e5)
#ax.set_xlim(3,50)
ax.semilogy()
ax.semilogx()
ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=6)
ax.legend([hhera,hthy], 
          [r'$\rm HERA$',r'$\rm JAM$'],
           fontsize = 10, frameon = 0)
ax.text(3.4,1.65e7, r'$x=2.8\cdot 10^{-5}\, (i=24)$', fontsize = 5)
ax.text(5.5,8.65e6, r'$x=4.6\cdot 10^{-5}$'         , fontsize = 5)
ax.text(8.0,4.6e6, r'$x=7.3\cdot 10^{-5}$'          , fontsize = 5)
ax.text(15.0,2.7e6, r'$x=1.2\cdot 10^{-4}$'         , fontsize = 5)
ax.text(22.0   ,1.4e6, r'$x=1.9\cdot 10^{-4}$'    , fontsize = 5)
ax.text(35.0   ,7.5e5, r'$x=3.1\cdot 10^{-4}$'    , fontsize = 5)
ax.text(45.0   ,3.9e5, r'$x=4.4\cdot 10^{-4}$'    , fontsize = 5)
ax.text(60.0   ,1.8e5, r'$x=6.2\cdot 10^{-4}$'    , fontsize = 5)
ax.text(80.0   ,8.5e4, r'$x=7.7\cdot 10^{-4}$'    , fontsize = 5)
ax.text(86.0   ,4e4  , r'$x=9.1 \cdot 10^{-4}$'   , fontsize = 5)
ax.text(160.0  ,2e4  , r'$x=0.0014$'              , fontsize = 5)
ax.text(260.0  ,1e4  , r'$x=0.0024$'              , fontsize = 5)
ax.text(440.0  ,4500 , r'$x=0.0035$'              , fontsize = 5)
ax.text(740.0  ,2000 , r'$x=0.0056$'              , fontsize = 5)
ax.text(1.5e3  ,970  , r'$x=0.0083$'              , fontsize = 5)
ax.text(1.7e3  ,410  , r'$x=0.013$'               , fontsize = 5)
ax.text(2.7e3  ,190  , r'$x=0.021$'               , fontsize = 5)
ax.text(3.9e3  ,77   , r'$x=0.032$'               , fontsize = 5)
ax.text(6.5e3  ,35   , r'$x=0.052$'               , fontsize = 5)
ax.text(1.15e4 ,19   , r'$x=0.084$'               , fontsize = 5)
ax.text(1.15e4 ,5.5  , r'$x=0.13$'                , fontsize = 5)
ax.text(1.15e4 ,2.25 , r'$x=0.18$'                , fontsize = 5)
ax.text(1.15e4 ,0.97 , r'$x=0.25$'                , fontsize = 5)
ax.text(1.15e4 ,0.25 , r'$x=0.40$'                , fontsize = 5)
ax.text(1.0e4 ,0.06 , r'$x=0.65 \, (i=0)$'       , fontsize = 5)

ax.set_xticks([1.0, 10.0, 1e2, 1e3, 1e4, 1e5,])
ax.set_xticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_tick_params(which = 'minor', length = 3)
ax.set_ylabel(r'$\sigma_r^{p,NC}$'+r'$\, \times\, 2^{\, i}$', size = 10)
ax.set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 10)

#ax.set_title(path)
py.tight_layout()
checkdir('gallery')
py.savefig('gallery/main04.pdf')






