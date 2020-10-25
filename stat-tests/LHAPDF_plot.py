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

lhafl = {'TBAR': -6, 'BBAR': -5, 'CBAR': -4, 'SBAR': -3, 'UBAR': -2, 'DBAR': -1,
            'GLUON': 0, 'D': 1, 'U': 2, 'S': 3, 'C': 4, 'B': 5, 'T': 6, 'PHT': 7}


Q=1
X1=np.logspace(-5,-1,101)[:100]
X2=np.linspace(0.1,1,100)
X=X1 #np.concatenate([X1,X2])
Nrep=101

PDFs_names = ["NNPDF31_nnlo_pch_as_0118",
        "NNPDF31_nnlo_pch_as_0118_rs_0.5",
        "NNPDF31_nnlo_pch_as_0118_rs_1.0"]

ref_PDF = 'NNPDF31_nnlo_pch_as_0118'

PDFs_colors = ["green",
               "blue",
               "red"]

PDFs_labels = [r'NNPDF3.1 (NNLO)',
               r'$(H_0):\,R_s=0.5$',
               r'$(H_1):\,R_s=1.0$']

fls = [lhafl['UBAR'],lhafl['DBAR'],lhafl['S'],lhafl['SBAR']]
fls_labels = [r'$\bar{u}$',r'$\bar{d}$',r'$s$',r'$\bar{s}$']
fig_composition= [2,2]

PDFs = {}
Rs = {}

for iPDF,PDF_name in enumerate(PDFs_names):
    PDFs[PDF_name]={}
    Rs[PDF_name] = []
    PDFs[PDF_name]['color']=PDFs_colors[iPDF]
    PDFs[PDF_name]['label'] = PDFs_labels[iPDF]
    for fl in range(-6,8):
        PDFs[PDF_name][fl]={'reps':[]}

    for irep in range(1,Nrep):
        LHAPDF = lhapdf.mkPDF(PDF_name, irep)

        for fl in range(-6, 8):
            rep = []
            for x in X:
                rep.append(LHAPDF.xfxQ(fl, x, Q))
                PDFs[PDF_name][fl]['reps'].append(rep)

for iPDF, PDF_name in enumerate(PDFs_names):
    Rs[PDF_name]=(np.array(PDFs[PDF_name][lhafl['SBAR']]['reps']) +
            np.array(PDFs[PDF_name][lhafl['S']]['reps']))/(np.array(PDFs[PDF_name][lhafl['UBAR']]['reps'])+np.array(PDFs[PDF_name][lhafl['DBAR']]['reps']))
        
Rs_median={}
Rs_low={}
Rs_up={}
for PDF_name in PDFs_names:
    for fl in range(-6,8):
        PDFs[PDF_name][fl]['low95cl'] = np.nanpercentile(PDFs[PDF_name][fl]['reps'], 5., axis=0)
        PDFs[PDF_name][fl]['up95cl'] = np.nanpercentile(PDFs[PDF_name][fl]['reps'], 95., axis=0)
        PDFs[PDF_name][fl]['median'] = np.median(PDFs[PDF_name][fl]['reps'], axis=0)

    Rs_median[PDF_name]=np.median(Rs[PDF_name], axis=0)
    Rs_low[PDF_name]=np.nanpercentile(Rs[PDF_name], 5., axis=0)
    Rs_up[PDF_name]=np.nanpercentile(Rs[PDF_name], 95., axis=0)


fig = gcf()
fig.suptitle(PDFs_names[0].replace("_","\_")+" at Q = "+str(Q)+" GeV", fontsize=14)

for ifl, fl in enumerate(fls):
    ax = py.subplot(fig_composition[0]*100+fig_composition[1]*10+ifl+1)
    for PDF_name in PDFs_names:

        ax.plot(X, PDFs[PDF_name][fl]['median'], ls='-', color=PDFs[PDF_name]['color'], lw=1.25,label=PDFs[PDF_name]['label'])

        ax.fill_between(X, PDFs[PDF_name][fl]['up95cl'], PDFs[PDF_name][fl]['low95cl'],
                        facecolor=PDFs[PDF_name]['color'], alpha=0.25, edgecolor=None, lw=1.25)
    
    ax.set_xscale('log')
    ax.set_xlabel('x')
    ax.set_ylabel('xf(x)')
    if ifl ==0:
        ax.legend(loc='best', title=fls_labels[ifl])
    else:
        props = dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.5)
        ax.text(0.9, 0.9, fls_labels[ifl], transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)

if not os.path.isdir("plots"):
    os.mkdir("plots")
                
resultpath = "plots/PDFs.pdf"
print resultpath+"... "
py.savefig(resultpath)
py.cla()
py.clf()

#------------- Rs
fig = gcf()
#fig.suptitle(PDFs_names[0].replace("_","\_")+r" $R_s = \frac{s+\bar{s}}{\bar{u}+\bar{d}}$ at Q = "+str(Q)+" GeV", fontsize=14)

np.seterr(divide='ignore', invalid='ignore')

#pdfbaseline = "NNPDF31_nnlo_pch_as_0118"
#title = scenario_choice+" Scenario"
#fig.suptitle(r'\hspace{-15pt}\textbf{'+pdfbaseline.replace("_", "\_") + r'}' , fontsize=19, y=0.993)

ax = py.subplot(111)
max_xbin = -12  # add [:max_xbin] to PDFs
for PDF_name in PDFs_names:

    #R_s = (np.array(PDFs[PDF_name][lhafl['S']]['median'])+np.array(PDFs[PDF_name][lhafl['SBAR']]['median']))/(np.array(PDFs[PDF_name][lhafl['UBAR']]['median'])+np.#array(PDFs[PDF_name][lhafl['DBAR']]['median']))       
#
    #R_s_up = (np.array(PDFs[PDF_name][lhafl['S']]['up95cl'])+np.array(PDFs[PDF_name][lhafl['SBAR']]['up95cl']))/(
    #    np.array(PDFs[PDF_name][lhafl['UBAR']]['up95cl'])+np.array(PDFs[PDF_name][lhafl['DBAR']]['up95cl']))
#
    #R_s_low = (np.array(PDFs[PDF_name][lhafl['S']]['low95cl'])+np.array(PDFs[PDF_name][lhafl['SBAR']]['low95cl']))/(
    #    np.array(PDFs[PDF_name][lhafl['UBAR']]['low95cl'])+np.array(PDFs[PDF_name][lhafl['DBAR']]['low95cl']))

    ax.plot(X, Rs_median[PDF_name], ls='-', color=PDFs[PDF_name]
            ['color'], lw=4.25, label=PDFs[PDF_name]['label'])

    ax.fill_between(X, Rs_up[PDF_name], Rs_low[PDF_name],
                    facecolor=PDFs[PDF_name]['color'], alpha=0.25, edgecolor=None, lw=1.25)

ax.set_xlim(1e-4,1e-1)
ax.set_ylim(0.,1.5)
ax.set_xscale('log')
ax.set_xlabel(r'{\boldmath $x$}', fontsize=45)
ax.xaxis.set_label_coords(0.9, +0.15)
#ax.set_ylabel(r'{\boldmath $R_s = \frac{s+\bar{s}}{\bar{u}+\bar{d}}$}',fontsize=25,rotation=0)
ax.set_ylabel(r'{\boldmath $R_s$}',fontsize=45,rotation=0)
ax.yaxis.set_label_coords(+0.1, 0.8)
#ax.legend(loc='best', title=r'NNPDF3.1 (nnlo)', fontsize=25)
#ax.legend(loc='best', fontsize=23)

#ax.text(0.45, 0.5, r'$Q='+str(Q)+"$ GeV, $90\%$ CL", transform=ax.transAxes, fontsize=20,
#        verticalalignment='top', bbox=props)

py.tight_layout()
         
resultpath = "plots/R_s.pdf"
print resultpath+"... "
py.savefig(resultpath)
py.cla()
py.clf()
#-------------
"""
fig = gcf()
fig.suptitle(PDFs_names[0].replace("_","\_")+r" $R_{f_v} = \frac{f_v}{f^{ref}_v}$ at Q = "+str(Q)+" GeV", fontsize=14)

np.seterr(divide='ignore', invalid='ignore')

max_xbin=-12
for PDF_name in PDFs_names:

    ax = py.subplot(211)
    R_uv = (np.array(PDFs[PDF_name][lhafl['U']]['median'][:max_xbin])-np.array(PDFs[PDF_name][lhafl['UBAR']]['median'][:max_xbin]))/(np.array(PDFs[ref_PDF][lhafl['U']]['median'][:max_xbin])-np.array(PDFs[ref_PDF][lhafl['UBAR']]['median'][:max_xbin])) 

    ax.plot(X[:max_xbin], R_uv, ls='-', color=PDFs[PDF_name]['color'], lw=1.25,label=PDFs[PDF_name]['label'])
    ax.set_xlabel('x')
    ax.set_ylabel(r'$R_{u_v}(x)$')
    ax.set_xscale('log')

    ax = py.subplot(212)
    R_dv = (np.array(PDFs[PDF_name][lhafl['D']]['median'][:max_xbin])-np.array(PDFs[PDF_name][lhafl['DBAR']]['median'][:max_xbin]))/(
        np.array(PDFs[ref_PDF][lhafl['D']]['median'][:max_xbin])-np.array(PDFs[ref_PDF][lhafl['DBAR']]['median'][:max_xbin]))

    ax.plot(X[:max_xbin], R_dv, ls='-', color=PDFs[PDF_name]['color'], lw=1.25,label=PDFs[PDF_name]['label'])
    ax.set_xlabel('x')
    ax.set_ylabel(r'$R_{d_v}(x)$')
    ax.set_xscale('log')


ax.legend(loc='best', title=r"$R_{f_v}$")
         
resultpath = "plots/R_v.pdf"
print resultpath+" saved... "
py.savefig(resultpath)
py.cla()
py.clf()
"""

print "DONE!!"

