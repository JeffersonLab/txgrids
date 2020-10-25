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
rc('xtick', labelsize=25)
rc('ytick', labelsize=25)

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

fls = [lhafl['F2NC']]
fls_labels = [r'$F_2^{NC}$']  # +90\%CL
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
#fig.suptitle(SFs_names[0].replace("_","\_")+" at Q = "+str(Q)+" GeV", fontsize=14)
pdfbaseline = "NNPDF31_nnlo_pch_as_0118"
fig.suptitle(r'\hspace{-15pt}\textbf{'+pdfbaseline.replace("_", "\_") + r'}' , fontsize=19, y=0.993)

for ifl, fl in enumerate(fls):
    ax = py.subplot(fig_composition[0]*100+fig_composition[1]*10+ifl+1)
    for SF_name in SFs_names:

        ax.plot(X, SFs[SF_name][fl]['median'], ls='-', color="blue", alpha=1.0, lw=2.25, label=r'$(H_0):\,Central$')
        ax.plot(X, np.ones(len(X))*100, ls='-', color="red", alpha=1.0, lw=2.25, label=r'$(H_1):\,Replica (k)$')

        #ax.fill_between(X, SFs[SF_name][fl]['up95cl'], SFs[SF_name][fl]['low95cl'],
        #                facecolor=SFs[SF_name]['color'], alpha=0.25, edgecolor=None, lw=2.25)
        
        for irep in range(len(SFs[SF_name][fl]['reps'])):
                ax.plot(X, SFs[SF_name][fl]['reps'][irep], ls='-', color="red", alpha=0.3, lw=0.5)
        
    
    ax.set_xscale('log')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(fls_labels[ifl])

    ax.set_xlim(1e-4,1)
    ax.set_ylim(0.0,0.6)

    ax.set_xscale('log')
    ax.set_xlabel(r'$x$',fontsize=25)
    ax.xaxis.set_label_coords(0.9, +0.1)
    ax.set_ylabel(fls_labels[ifl], fontsize=25, rotation=0)
    ax.yaxis.set_label_coords(+0.1, 0.8)

    ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
    ax.set_xticklabels([r'', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    #ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    #ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)


    if ifl ==0:
        ax.legend(loc='lower left', fontsize=25)
    else:
        props = dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.5)
        ax.text(0.9, 0.9, fls_labels[ifl], transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

if not os.path.isdir("plots"):
    os.mkdir("plots")
                
resultpath = "plots/F2.pdf"
print resultpath+"... "
py.tight_layout()
py.savefig(resultpath)
py.cla()
py.clf()

print "DONE!!"

