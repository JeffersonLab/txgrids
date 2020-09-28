'''
 # @ Author: Rabah Abdul Khalek
 # @ Create Time: 2020-04-25 16:33:25
 # @ Modified by: Rabah Abdul Khalek
 # @ Modified time: 2020-04-25 19:20:16
 # @ Description: New ideas for impact studies at  EIC
 '''

import os
import sys
import numpy as np
import matplotlib.pyplot as py
from pylab import *
from matplotlib import rc
import collections
from scipy.integrate import quad
import yaml
import lhapdf
cfg = lhapdf.getConfig()
cfg.set_entry("Verbosity", 0)
sys.path.append(os.getcwd()+"/../")
rc('font', **{'family': 'sans-serif', 'sans-serif': []})
rc('text', usetex=True)

lhafl = {'F2NC': 908, 'FLNC': 909, 'F3NC':910}


Q=1
X1=np.logspace(-4,-1,101)[:100]
X2=np.linspace(0.1,1,100)
X=np.concatenate([X1,X2])
Nrep=101

SFs_names = ["NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF"]
            #,"NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF_FLNC_90CL_LOW", "NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF_FLNC_90CL_UP"]

ref_SF = 'NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'

SFs_colors = ["green","blue","red"]

SFs_labels = [r'Reference',"low","up"]

fls = [lhafl['FLNC']]
fls_labels = [r'$F_L^{NC}+90\%CL$']
fig_composition= [1,1]

SFs = {}

for iSF,SF_name in enumerate(SFs_names):
    SFs[SF_name]={}
    SFs[SF_name]['color']=SFs_colors[iSF]
    SFs[SF_name]['label'] = SFs_labels[iSF]
    for fl in fls:
        SFs[SF_name][fl]={'reps':[]}

    for irep in range(1,Nrep):
        LHAPDF = lhapdf.mkPDF(SF_name, irep)

        for fl in fls:
            rep = []
            for x in X:
                rep.append(LHAPDF.xfxQ(fl, x, Q))
                SFs[SF_name][fl]['reps'].append(rep)

for SF_name in SFs_names:
    for fl in fls:
        SFs[SF_name][fl]['low95cl'] = np.nanpercentile(SFs[SF_name][fl]['reps'], 5., axis=0)
        SFs[SF_name][fl]['up95cl'] = np.nanpercentile(SFs[SF_name][fl]['reps'], 95., axis=0)
        SFs[SF_name][fl]['median'] = np.median(SFs[SF_name][fl]['reps'], axis=0)

fig = gcf()
fig.suptitle(SFs_names[0].replace("_","\_")+" at Q = "+str(Q)+" GeV", fontsize=14)

for ifl, fl in enumerate(fls):
    ax = py.subplot(fig_composition[0]*100+fig_composition[1]*10+ifl+1)
    for SF_name in SFs_names:

        ax.plot(X, SFs[SF_name][fl]['median'], ls='-', color=SFs[SF_name]['color'], lw=1.25,label=SFs[SF_name]['label'])

        ax.fill_between(X, SFs[SF_name][fl]['up95cl'], SFs[SF_name][fl]['low95cl'],
                        facecolor=SFs[SF_name]['color'], alpha=0.25, edgecolor=None, lw=1.25)
    
    ax.set_xscale('log')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(fls_labels[ifl])
    if ifl ==0:
        ax.legend(loc='best', title=fls_labels[ifl])
    else:
        props = dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.5)
        ax.text(0.9, 0.9, fls_labels[ifl], transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

if not os.path.isdir("plots"):
    os.mkdir("plots")
                
resultpath = "plots/FL.pdf"
print resultpath+"... "
py.savefig(resultpath)
py.cla()
py.clf()

print "DONE!!"

