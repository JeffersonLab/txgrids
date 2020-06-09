'''
 # @ Author: Rabah Abdul Khalek
 # @ Create Time: 2020-04-24 14:02:30
 # @ Modified by: Rabah Abdul Khalek
 # @ Modified time: 2020-04-24 15:04:17
 # @ Description: New ideas for impact studies at  EIC
 '''

import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory.tools import save, load, convert_lum
from theory.mceg  import MCEG
from ndtest import ndtest

#--matplotlib
import matplotlib
import matplotlib.pyplot as plt
import pylab  as py
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

from scipy import stats
import pandas as pd
from collections import OrderedDict, Counter

from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

def units(units):
    one=0.3893793721 #--GeV2 mbarn
    if   units=='GeV^-2':return 1
    elif units=='fb'    :return one*1e12 
    else: sys.exit('%s conversion not available')

wdir='.stats-tests'

if len(sys.argv) < 1:
    print "usage: ./stat-tests.py <Luminosity (fb^-1)>"
    exit(1)

lum_arg = float(sys.argv[1])

#to tweak:
#-----------------------------------------------

#MCEG integration accuracy
neval = 10000  # (neval/2) ~ batch_size
nitn = 1000
min_Q = 0.1

#SFs
min_SF = 'JAM4EIC'
max_SF = 'JAM4EIC'

#min_SF = 'NNPDF31_nnlo_pch_as_0118_rs_0.5_SF'
#max_SF = 'NNPDF31_nnlo_pch_as_0118_rs_1.0_SF'

seeds = {'min': 1, 'max': 2}

#---binning
NQ2bins = 20
print "NQ2bins = ", NQ2bins
NXbins = 20
print "NXbins = ", NXbins
choice_bins = 'min'
#---

#--physical params and cuts
rs= 140.7
lum=str(lum_arg)+':fb-1'
lum_label = r'$\mathcal{L} = '+str(lum_arg)+r'\,fb^{-1}$'
sign=1 #--electron=1 positron=-1
cuts_label = r"$ W^2 > 10\,\,,\,\, Q^2> 1$"
def veto00(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1
#-----------------------------------------------

#--lhapdf set and stf idx
s = {'min':{}, 'max':{}}

s['min']['tabname'] = min_SF
s['min']['iset'] = 0
s['min']['iF2'] = 908
s['min']['iFL'] = 909
s['min']['iF3'] = 910

s['max']['tabname'] = max_SF
s['max']['iset'] = 0
s['max']['iF2'] = 908
s['max']['iFL'] = 909
s['max']['iF3'] = 910

#common keys
for key in s.keys():
    s[key]['wdir']      = wdir
    s[key]['sign']      = sign
    s[key]['rs'] = rs
    s[key]['fname'] = 'mceg00'
    s[key]['veto'] = veto00
    s[key]['fdata'] = wdir+"/"+s[key]['tabname']+"_data_lum"+str(lum_arg)+"_seed"+str(seeds[key])+".po"

#--generate events
#output keys: ['Y', 'X', 'Q2', 'W', 'rs']

for key in s.keys():
    np.random.seed(seeds[key])
    if not os.path.isfile(s[key]['fdata']):
        mceg = MCEG(**(s[key]))
        mceg.buil_mceg(neval=neval, nitn=nitn, min_Q=min_Q)
        s[key]['tot_xsec'] = mceg.get_xsectot()
        s[key]['tot_xsec'] *= units('fb')

        s[key]['var_xsec'] = mceg.get_xsecvar()
        s[key]['var_xsec'] = s[key]['var_xsec']**0.5*units('fb')

        s[key]['quality_xsec'] = mceg.get_quality()

        Nevents = int(lum_arg*s[key]['tot_xsec'])
        
        s[key].update(mceg.gen_events(Nevents))

        save(s[key], s[key]['fdata'])
    else:
        s[key] = load(s[key]['fdata'])

np.set_printoptions(threshold=sys.maxsize)

keys_of_interest = ['W','X','Q2']

for key1 in s.keys():
    for key2 in keys_of_interest:
        s[key1][key2]=np.array(s[key1][key2])
#------
#--------------------------------------------------------------------------------------------------------------------


Q2bins_smin = np.logspace(np.log10(np.min(s['min']['Q2'])), np.log10(np.max(s['min']['Q2'])),NQ2bins+1)
Q2bins_smax = np.logspace(np.log10(np.min(s['max']['Q2'])), np.log10(np.max(s['max']['Q2'])),NQ2bins+1)

Xbins_smin = np.logspace(np.log10(np.min(s['min']['X'])), np.log10(np.max(s['min']['X'])),NXbins+1)
Xbins_smax = np.logspace(np.log10(np.min(s['max']['X'])), np.log10(np.max(s['max']['X'])),NXbins+1)

if choice_bins == 'min':
    Q2bins = Q2bins_smin
    Xbins = Xbins_smin
    resultpath = 'plots/chi2-test_perxQ2_'+str(NXbins)+"-"+str(NQ2bins)+'sminbins_'+str(Nevents/1000)+'k-events.pdf'
elif choice_bins == 'max':
    Q2bins = Q2bins_smax
    Xbins = Xbins_smax
    resultpath = 'plots/chi2-test_perxQ2_'+str(NXbins)+"-"+str(NQ2bins)+'smaxbins_'+str(Nevents/1000)+'k-events.pdf'


non_empty = {'min': np.zeros(NXbins*NQ2bins)>1, 'max': np.zeros(NXbins*NQ2bins)>1}
covmat_xsecs = {'min': np.zeros((NXbins*NQ2bins, NXbins*NQ2bins)), 'max': np.zeros((NXbins*NQ2bins, NXbins*NQ2bins))}

hist_weights = {'min': np.zeros(NXbins*NQ2bins), 'max': np.zeros(NXbins*NQ2bins)}
hist_xsecs = {'min': np.zeros(NXbins*NQ2bins), 'max': np.zeros(NXbins*NQ2bins)}
hist_N = {'min': np.zeros(NXbins*NQ2bins), 'max': np.zeros(NXbins*NQ2bins)}

#to get:
#Nevents: hist_weights['min']*s['max']['tot_xsec']*lum (per bin)
#tot xsec: s['min']['tot_xsec'] (in fb)
#normalized xsec: hist_weights['min'

for key in hist_weights.keys():
    for iQ2 in range(0, NQ2bins):
        for iX in range(0, NXbins):

            x_mask  = np.where((s[key]['X'] > Xbins[iX]) & (s[key]['X'] < Xbins[iX+1]),True,False)
            Q2_mask = np.where((s[key]['Q2'] > Q2bins[iQ2]) & (s[key]['Q2'] < Q2bins[iQ2+1]),True,False)
            mask=x_mask*Q2_mask
            
            weight = np.sum(s[key]['W'][mask])
            xsec = weight*s[key]['tot_xsec'] # [fb]
            N = float(int(xsec*lum_arg))

            hist_weights[key][iQ2*NXbins+iX] = weight
            hist_xsecs[key][iQ2*NXbins+iX] = xsec
            hist_N[key][iQ2*NXbins+iX] = N

            if N==0:
                non_empty[key][iQ2*NXbins+iX] = False
                covmat_xsecs[key][iQ2*NXbins+iX][iQ2*NXbins+iX] = 0.
            else:
                non_empty[key][iQ2*NXbins+iX] = True
                covmat_xsecs[key][iQ2*NXbins+iX][iQ2*NXbins +
                                                 iX] = (np.sqrt(N)/lum_arg)**2 + (weight*s[key]['var_xsec'])**2
                #covmat_weights[key][iQ2*NXbins+iX][iQ2*NXbins+iX] = (np.sqrt(N)/(s[key]['tot_xsec']*lum))**2
            """
            print key
            print 'xsec = ', xsec
            print 'Lum =',lum
            print 'Nevents = ',N
            print 'delta_xsec = ', (np.sqrt(N)/lum)
            print ' '
            """

tot_non_empty = non_empty['min']+non_empty['max']

map_non_empty={}
count=0
for iQ2 in range(0, NQ2bins):
    for iX in range(0, NXbins):
        if tot_non_empty[iQ2*NXbins+iX]:
            map_non_empty[iQ2*NXbins+iX]=count
            count+=1

N_bins = np.count_nonzero(tot_non_empty)

tot_covmat_xsecs = np.matrix(np.zeros((N_bins, N_bins)))

weights = {}
xsecs = {}
for key in hist_xsecs.keys():
    xsecs[key] = np.matrix(hist_xsecs[key][tot_non_empty])
    covmat_xsecs[key] = covmat_xsecs[key][tot_non_empty, :]
    covmat_xsecs[key] = covmat_xsecs[key][:, tot_non_empty]

    covmat_xsecs[key] = np.matrix(covmat_xsecs[key])
    tot_covmat_xsecs += covmat_xsecs[key]

inv_covmat_xsecs = np.linalg.inv(tot_covmat_xsecs)

hist_chi2s = np.zeros(NXbins*NQ2bins)
hist_Xs = np.zeros(NXbins*NQ2bins)
hist_Q2s = np.zeros(NXbins*NQ2bins)
for iQ2 in range(0, NQ2bins):
    for iX in range(0, NXbins):
        if tot_non_empty[iQ2*NXbins+iX]:
            ind = map_non_empty[iQ2*NXbins+iX]

            max_d = xsecs['max'][0,ind]
            min_d = xsecs['min'][0,ind]

            hist_chi2s[iQ2*NXbins+iX] = (min_d-max_d)*inv_covmat_xsecs[ind,ind]*(min_d-max_d)
        
        hist_Xs[iQ2*NXbins+iX] = Xbins[iX]+(Xbins[iX+1]-Xbins[iX])/2.
        hist_Q2s[iQ2*NXbins+iX] = Q2bins[iQ2]+(Q2bins[iQ2+1]-Q2bins[iQ2])/2.


tot_chi2 = (xsecs['min']-xsecs['max'])*inv_covmat_xsecs*(xsecs['min']-xsecs['max']).T
tot_chi2 /= N_bins

nrows, ncols = 1, 3
gs = gridspec.GridSpec(nrows, ncols)
fig = py.figure(figsize=(ncols*5, nrows*3.5))

fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+s['min']['tabname'].replace("_", "\_") + r'}$,' +
             r' $\textrm{max: '+s['max']['tabname'].replace("_", "\_") + r'}$, '+
             lum_label, fontsize=10, y=0.98)

fig.subplots_adjust(top=0.85, bottom=0.15)

ax = py.subplot(gs[0])

#---- chi2 hist
h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_chi2s, bins=[np.log(Xbins), np.log(Q2bins)], cmap='hot_r')

ax.title.set_text(r'$\chi^2_{tot}/N_{bins}$ = '+' %0.4e'%tot_chi2)
cbar = plt.colorbar()
cbar.ax.set_xlabel(r"$\chi^2_{bin}$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$Q^2$")

ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])

#----- min xsec hist
ax = py.subplot(gs[1])
ax.title.set_text(r'$\sigma^{min}_{tot}$'+' = %0.2e' % (s['min']['tot_xsec']) + r'$\pm$'+'%0.2e' % (s['min']['var_xsec'])+r' [fb]')
h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_xsecs['min'], bins=[np.log(Xbins), np.log(Q2bins)], cmap='hot_r')
cbar = plt.colorbar()
cbar.ax.set_xlabel(r"$\frac{d\sigma_{min}}{dxdQ^2}$")
#cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{min}}{dxdQ^2}$")
ax.set_xlabel(r"$x$")

ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])

#----- max Nevents hist
ax = py.subplot(gs[2])
ax.title.set_text(r'$\sigma^{max}_{tot}$'+' = %0.2e' % (s['max']['tot_xsec']) + r'$\pm$'+'%0.2e' % (s['max']['var_xsec'])+r' [fb]')
h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_xsecs['max'], bins=[np.log(Xbins), np.log(Q2bins)], cmap='hot_r')
cbar = plt.colorbar()
cbar.ax.set_xlabel(r"$\frac{d\sigma_{max}}{dxdQ^2}$")
#cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
ax.set_xlabel(r"$x$")

ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])

plt.savefig("chi2-test-"+str(lum_arg)+"fb-1.pdf")




