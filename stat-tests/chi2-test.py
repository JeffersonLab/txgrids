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
from expdata.driver import gen_grid

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

def GetEvents(lum_arg,veto,sign):
    """
    sign: --electron=1 positron=-1

    """
    Sample = {'min':{}, 'max':{}}

    Sample['min']['tabname'] = min_SF
    Sample['min']['iset'] = 0
    Sample['min']['iF2'] = 908
    Sample['min']['iFL'] = 909
    Sample['min']['iF3'] = 910

    Sample['max']['tabname'] = max_SF
    Sample['max']['iset'] = 0
    Sample['max']['iF2'] = 908
    Sample['max']['iFL'] = 909
    Sample['max']['iF3'] = 910

    #common keys
    for key in Sample.keys():
        Sample[key]['wdir']      = wdir
        Sample[key]['sign']      = sign
        Sample[key]['rs'] = rs
        Sample[key]['fname'] = 'mceg00'
        Sample[key]['veto'] = veto
        Sample[key]['fdata'] = wdir+"/"+Sample[key]['tabname']+"_data_lum"+str(lum_arg)+"_seed"+str(seeds[key])+".po"

    #--generate events
    #output keys: ['Y', 'X', 'Q2', 'W', 'rs']
    print("--Getting Events--")
    for key in Sample.keys():
        np.random.seed(seeds[key])
        if not os.path.isfile(Sample[key]['fdata']):
            mceg = MCEG(**(Sample[key]))
            mceg.buil_mceg(neval=neval, nitn=nitn, min_Q=min_Q)

            Sample[key]['tot_xsec'] = mceg.get_xsectot()
            Sample[key]['tot_xsec'] *= units('fb')

            Sample[key]['var_xsec'] = mceg.get_xsecvar()
            Sample[key]['var_xsec'] = Sample[key]['var_xsec']**0.5*units('fb')

            Sample[key]['quality_xsec'] = mceg.get_quality()

            Nevents = int(lum_arg*Sample[key]['tot_xsec'])
            
            Sample[key]['Nevents'] = Nevents
            
            Sample[key].update(mceg.gen_events(Nevents))

            save(Sample[key], Sample[key]['fdata'])
        else:
            Sample[key] = load(Sample[key]['fdata'])

    np.set_printoptions(threshold=sys.maxsize)

    keys_of_interest = ['W','X','Q2']

    for key1 in Sample.keys():
        for key2 in keys_of_interest:
            Sample[key1][key2]=np.array(Sample[key1][key2])

    return Sample

def GetBins():
    Bins = {}
    kin0, kin = gen_grid("../expdata/")
    Bins['x'] = sorted(set(np.array(kin)[:, 0]))
    Bins['Q2'] = sorted(set(np.array(kin)[:, 1]))

    Bins['x_cv'] = sorted(set(np.array(kin0)[:, 0]))[:-1]
    Bins['Q2_cv'] = sorted(set(np.array(kin0)[:, 1]))[:-1]

    Bins['NQ2'] = len(Bins['Q2'])-1
    Bins['Nx'] = len(Bins['x'])-1

    """ #user defined binning
    Bins['NQ2'] = 50
    print "Bins['NQ2'] = ", Bins['NQ2']
    Bins['Nx'] = 50
    print "Bins['Nx'] = ", Bins['Nx']
    choice_bins = 'min'

    Q2bins_smin = np.logspace(np.log10(np.min(Sample['min']['Q2'])), np.log10(np.max(Sample['min']['Q2'])),Bins['NQ2']+1)
    Q2bins_smax = np.logspace(np.log10(np.min(Sample['max']['Q2'])), np.log10(np.max(Sample['max']['Q2'])),Bins['NQ2']+1)

    Xbins_smin = np.logspace(np.log10(np.min(Sample['min']['X'])), np.log10(np.max(Sample['min']['X'])),Bins['Nx']+1)
    Xbins_smax = np.logspace(np.log10(np.min(Sample['max']['X'])), np.log10(np.max(Sample['max']['X'])),Bins['Nx']+1)

    if choice_bins == 'min':
        Bins['Q2'] = Q2bins_smin
        Bins['x'] = Xbins_smin
    elif choice_bins == 'max':
        Bins['Q2'] = Q2bins_smax
        Bins['x'] = Xbins_smax
    """
    return Bins

def GetCrossSections(Bins,sys_path="../expdata/data/xQ2binTable-xiaoxuan-060220+syst.npy"):

    #--- systematics uncertainties
    sys_bins = np.load(sys_path)

    #initialization
    hist_CrossSection={}
    hist_CrossSection['non_empty'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']) > 1, 'max': np.zeros(Bins['Nx']*Bins['NQ2']) > 1}
    hist_CrossSection['weights'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}
    hist_CrossSection['xsecs'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}
    hist_CrossSection['Nevents'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}
    hist_CrossSection['statunc'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}
    hist_CrossSection['MCunc'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}
    hist_CrossSection['sysunc'] = {'min': np.zeros(Bins['Nx']*Bins['NQ2']), 'max': np.zeros(Bins['Nx']*Bins['NQ2'])}

    #to get:
    #Sample[key]['Nevents']: hist_CrossSection['weights']['min']*Sample['max']['tot_xsec']*lum (per bin)
    #tot xsec: Sample['min']['tot_xsec'] (in fb)
    #normalized xsec: hist_CrossSection['weights']['min']

    print("--Filling hist_CrossSection--")
    for key in hist_CrossSection['weights'].keys():
        for iQ2 in range(0, Bins['NQ2']):
            for iX in range(0, Bins['Nx']):

                x_mask  = np.where((Sample[key]['X'] > Bins['x'][iX]) & (Sample[key]['X'] < Bins['x'][iX+1]),True,False)
                Q2_mask = np.where((Sample[key]['Q2'] > Bins['Q2'][iQ2]) & (Sample[key]['Q2'] < Bins['Q2'][iQ2+1]),True,False)
                mask=x_mask*Q2_mask
                weight = np.sum(Sample[key]['W'][mask])
                xsec = weight*Sample[key]['tot_xsec'] # [fb]
                N = float(int(xsec*lum_arg))
                hist_CrossSection['weights'][key][iQ2*Bins['Nx']+iX] = weight
                hist_CrossSection['xsecs'][key][iQ2*Bins['Nx']+iX] = xsec
                hist_CrossSection['Nevents'][key][iQ2*Bins['Nx']+iX] = N
                if N==0:
                    hist_CrossSection['non_empty'][key][iQ2*Bins['Nx']+iX] = False
                    hist_CrossSection['statunc'][key][iQ2*Bins['Nx']+iX] = 0.
                    hist_CrossSection['MCunc'][key][iQ2*Bins['Nx']+iX] = 0.
                    hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = 0.
                else:
                    hist_CrossSection['non_empty'][key][iQ2*Bins['Nx']+iX] = True
                    hist_CrossSection['statunc'][key][iQ2*Bins['Nx']+iX] = np.sqrt(N)/lum_arg
                    hist_CrossSection['MCunc'][key][iQ2*Bins['Nx']+iX] = weight*Sample[key]['var_xsec']

                    relative_sys = np.where(
                        (sys_bins[:, 0] == Bins['x_cv'][iX]) & (sys_bins[:, 1] == Bins['Q2_cv'][iQ2]))
                    if np.any(relative_sys): #if bin has sys uncertainty
                        relative_sys = sys_bins[relative_sys][0, 2]
                    else: #if bin doesn't have sys uncertainty
                        relative_sys = 0

                    hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = xsec*relative_sys

    hist_CrossSection['non_empty']['total'] = hist_CrossSection['non_empty']['min']+hist_CrossSection['non_empty']['max']

    hist_CrossSection['non_empty']['map'] = {}
    count = 0
    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if hist_CrossSection['non_empty']['total'][iQ2*Bins['Nx']+iX]:
                hist_CrossSection['non_empty']['map'][iQ2 *
                                                      Bins['Nx']+iX] = count
                count += 1

    Bins['total'] = np.count_nonzero(hist_CrossSection['non_empty']['total'])

    return hist_CrossSection

def GetCovMat(Bins, hist_CrossSection):

    CovMat = {'xsecs': {'min': np.zeros((Bins['Nx']*Bins['NQ2'], Bins['Nx']*Bins['NQ2'])), 'max': np.zeros(
        (Bins['Nx']*Bins['NQ2'], Bins['Nx']*Bins['NQ2'])), 'total': np.zeros((Bins['Nx']*Bins['NQ2'], Bins['Nx']*Bins['NQ2']))}}
        
    print("--Building covariance matrix--")
    for key in hist_CrossSection['weights'].keys():
        for iQ2 in range(0, Bins['NQ2']):
            for iX in range(0, Bins['Nx']):

                if hist_CrossSection['Nevents'][key][iQ2*Bins['Nx']+iX] == 0:
                    CovMat['xsecs'][key][iQ2*Bins['Nx']+iX][iQ2*Bins['Nx']+iX] = 0.
                else:
                    CovMat['xsecs'][key][iQ2*Bins['Nx']+iX][iQ2*Bins['Nx']+iX] = (hist_CrossSection['statunc'][key][iQ2*Bins['Nx']+iX])**2 + (
                        hist_CrossSection['MCunc'][key][iQ2*Bins['Nx']+iX])**2 + (hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX])**2
                    #covmat_weights[key][iQ2*Bins['Nx']+iX][iQ2*Bins['Nx']+iX] = (np.sqrt(N)/(Sample[key]['tot_xsec']*lum))**2

    CovMat['xsecs']['total'] = np.matrix(np.zeros((Bins['total'], Bins['total'])))

    for key in hist_CrossSection['xsecs'].keys():
        CovMat['xsecs'][key] = CovMat['xsecs'][key][hist_CrossSection['non_empty']['total'], :]
        CovMat['xsecs'][key] = CovMat['xsecs'][key][:, hist_CrossSection['non_empty']['total']]

        CovMat['xsecs'][key] = np.matrix(CovMat['xsecs'][key])
        CovMat['xsecs']['total'] += CovMat['xsecs'][key]
    
    return CovMat

def units(units):
    one=0.3893793721 #--GeV2 mbarn
    if   units=='GeV^-2':return 1
    elif units=='fb'    :return one*1e12 
    else: sys.exit('%s conversion not available')


if __name__ == "__main__":
        
    wdir='.stats-tests'

    if len(sys.argv) < 1:
        print "usage: ./stat-tests.py <Luminosity (fb^-1)>"
        exit(1)

    lum_arg = float(sys.argv[1])

    #--physical params and cuts
    #lum = str(lum_arg)+':fb-1'
    lum_label = r'$\mathcal{L} = '+str(lum_arg)+r'\,fb^{-1}$'
    cuts_label = r"$ W^2 > 10\,\,,\,\, Q^2> 1$"

    #to tweak:
    #-----------------------------------------------

    #MCEG integration accuracy
    neval = 10000  # default should be: 10000  # (neval/2) ~ batch_size
    nitn = 1000  # default should be: 1000
    min_Q = 0.1

    #SFs
    min_SF = 'NNPDF31_nnlo_pch_as_0118_rs_0.5_SF'
    max_SF = 'NNPDF31_nnlo_pch_as_0118_rs_1.0_SF'

    seeds = {'min': 3, 'max': 4}

    rs = 140.7
    #-----------------------------------------------

    #--- acceptance/purity based binning
    Bins=GetBins()

    #--- kinematic cuts
    def veto00(x,y,Q2,W2):
        if   W2 < 10       : return 0
        elif Q2 < 1        : return 0
        else               : return 1
        
    #--Getting events
    Sample = GetEvents(lum_arg,veto00,sign=1)

    #--Getting CrossSections
    hist_CrossSection = GetCrossSections(Bins,sys_path="../expdata/data/xQ2binTable-xiaoxuan-060220+syst.npy")

    #--Getting CovMat
    CovMat = GetCovMat(Bins, hist_CrossSection)
    inv_covmat_xsecs = np.linalg.inv(CovMat['xsecs']['total'])

    hist_chi2s = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_zscores = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_Xs = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_Q2s = np.zeros(Bins['Nx']*Bins['NQ2'])

    xsecs = {}
    for key in hist_CrossSection['xsecs'].keys():
        xsecs[key] = np.matrix(hist_CrossSection['xsecs'][key][hist_CrossSection['non_empty']['total']])

    print("--Filling histograms--")
    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if hist_CrossSection['non_empty']['total'][iQ2*Bins['Nx']+iX]:
                ind = hist_CrossSection['non_empty']['map'][iQ2*Bins['Nx']+iX]

                max_d = xsecs['max'][0,ind]
                min_d = xsecs['min'][0,ind]

                hist_chi2s[iQ2*Bins['Nx']+iX] = (min_d-max_d)*inv_covmat_xsecs[ind,ind]*(min_d-max_d)
                hist_zscores[iQ2*Bins['Nx']+iX] = np.sqrt(hist_chi2s[iQ2*Bins['Nx']+iX])
                if hist_zscores[iQ2*Bins['Nx']+iX]<1:
                    hist_zscores[iQ2*Bins['Nx']+iX]=-1
                else:
                    hist_zscores[iQ2*Bins['Nx']+iX] = int(hist_zscores[iQ2*Bins['Nx']+iX])

            hist_Xs[iQ2*Bins['Nx']+iX] = Bins['x_cv'][iX]
            hist_Q2s[iQ2*Bins['Nx']+iX] = Bins['Q2_cv'][iQ2]

            ##hist_Xs[iQ2*Bins['Nx']+iX] = Bins['x'][iX]+(Bins['x'][iX+1]-Bins['x'][iX])/2.
            ##hist_Q2s[iQ2*Bins['Nx']+iX] = Bins['Q2'][iQ2]+(Bins['Q2'][iQ2+1]-Bins['Q2'][iQ2])/2.

    tot_chi2 = (xsecs['min']-xsecs['max'])*inv_covmat_xsecs*(xsecs['min']-xsecs['max']).T
    tot_chi2 /= Bins['total']

    hist_zscores[hist_zscores >= 5] = 5

    #replace all 0 bins with nan for the hist to be properly colored:
    for key in hist_CrossSection['xsecs'].keys():
        hist_CrossSection['weights'][key][hist_CrossSection['weights'][key]==0] = np.nan
        hist_CrossSection['xsecs'][key][hist_CrossSection['xsecs'][key]==0] = np.nan
        hist_chi2s[hist_chi2s==0] = np.nan
        hist_zscores[hist_zscores == 0] = np.nan

    hist_zscores[hist_zscores == -1] = 0

    title_xpos = 0.5
    title_ypos = 1.04
    nrows, ncols = 2, 4
    widths=[2,1,1,1]
    heights=[1,1]
    gs = gridspec.GridSpec(nrows, ncols, width_ratios=widths,
                        height_ratios=heights)
    fig = py.figure(figsize=(ncols*5, nrows*3.5))#, constrained_layout=True)

    fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+Sample['min']['tabname'].replace("_", "\_") + r'}$,' +
                r' $\textrm{max: '+Sample['max']['tabname'].replace("_", "\_") + r'}$, '+
                lum_label, fontsize=10, y=0.98)

    fig.subplots_adjust(top=0.85, bottom=0.15)

    ax = py.subplot(gs[:,0])

    print("--Plotting--")
    print("--hist1--")

    #---- chi2 hist
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_chi2s, bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis',norm=matplotlib.colors.LogNorm())

    ax.text(title_xpos, title_ypos+0.02, r'$\chi^2_{tot}/N_{bins}$ = '+' %0.4e' % tot_chi2,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
        fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\chi^2_{bin}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in",length=7)
    ax.tick_params(which='minor', direction="in",length=4)

    #----- max xsec hist
    print("--hist2--")
    ax = py.subplot(gs[0,1])
    ax.text(title_xpos, title_ypos, r'$\sigma^{max}_{tot}$'+' = %0.2e' % (Sample['max']['tot_xsec']) + r'$\pm$'+'%0.2e' % (Sample['max']['var_xsec'])+r' [fb]',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)

    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['xsecs']['max'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())

    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\frac{d\sigma_{max}}{dxdQ^2}$")

    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    #ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    #----- min xsec hist
    print("--hist3--")
    ax = py.subplot(gs[1,1])
    ax.text(title_xpos, title_ypos, r'$\sigma^{min}_{tot}$'+' = %0.2e' % (Sample['min']['tot_xsec']) + r'$\pm$'+'%0.2e' % (Sample['min']['var_xsec'])+r' [fb]',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['xsecs']['min'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\frac{d\sigma_{min}}{dxdQ^2}$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{min}}{dxdQ^2}$")
    ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    #----- max sys unc hist
    print("--hist4--")
    ax = py.subplot(gs[0,2])
    ax.text(title_xpos, title_ypos, r'Sys Unc (max)', #r'$\delta^{max,sys}_{tot}$'+' = %0.2e' %
            #((np.sqrt(Sample['max']['Nevents'])/lum_arg)
            # * 100./Sample['max']['tot_xsec']) + r' \%',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['sysunc']['max']*100./hist_CrossSection['xsecs']['max'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{max}_{sys} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    #ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    #----- min sys unc hist
    print("--hist5--")
    ax = py.subplot(gs[1,2])
    ax.text(title_xpos, title_ypos, r'Sys Unc (min)',#r'$\delta^{min,sys}_{tot}$'+' = %0.2e' % ((np.sqrt(Sample["min"]['Nevents'])/lum_arg) * 100./Sample['min']['tot_xsec']) + r' \%',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['sysunc']['min']*100./hist_CrossSection['xsecs']['min'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{min}_{sys} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    #----- max stat unc hist
    print("--hist6--")
    ax = py.subplot(gs[0,3])
    ax.text(title_xpos, title_ypos, r'Stat Unc (max)', #r'$\delta^{max,stat}_{tot}$'+' = %0.2e' %
            #((np.sqrt(Sample['max']['Nevents'])/lum_arg)
            # * 100./Sample['max']['tot_xsec']) + r' \%',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['statunc']['max']*100./hist_CrossSection['xsecs']['max'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{max}_{stat} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    #ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    #----- min stat unc hist
    print("--hist7--")
    ax = py.subplot(gs[1,3])
    ax.text(title_xpos, title_ypos, r'Stat Unc (min)', #r'$\delta^{min,stat}_{tot}$'+' = %0.2e' %
            #((np.sqrt(Sample["min"]['Nevents'])/lum_arg)
            # * 100./Sample['min']['tot_xsec']) + r' \%',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['statunc']['min']*100./hist_CrossSection['xsecs']['min'], bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{min}_{stat} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)
    """
    #----- max MC unc hist
    ax = py.subplot(gs[0,3])
    ax.set_title(r'$\delta^{max,MC}_{tot}$'+' = %0.2e' %
                    ((Sample['max']['var_xsec'])*100./Sample['max']['tot_xsec']) + r' \%')
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['MCunc']['max']*100./hist_CrossSection['xsecs']['max'], bins=[np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{max}_{MC} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])


    #----- min MC unc hist
    ax = py.subplot(gs[1,3])
    ax.set_title(r'$\delta^{min,MC}_{tot}$'+' = %0.2e' %
                    ((Sample['min']['var_xsec'])*100./Sample['min']['tot_xsec']) + r' \%')
    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_CrossSection['MCunc']['min']*100./hist_CrossSection['xsecs']['min'], bins=[np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{min}_{MC} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])

    """

    plt.savefig("chi2-test-"+str(lum_arg)+"fb-1.pdf")
    plt.cla()
    plt.clf()


    #---- Chi2 hist alone with values written per bin
    fig, ax = plt.subplots()
    ax.set_aspect("equal")


    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_chi2s, bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\chi^2_{bin}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$")
    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if hist_CrossSection['non_empty']['total'][iQ2*Bins['Nx']+iX]:
                ax.text(np.log(hist_Xs[iQ2*Bins['Nx']+iX]), np.log(hist_Q2s[iQ2*Bins['Nx']+iX]), '{:.2f}'.format(hist_chi2s[iQ2*Bins['Nx']+iX]),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=6)

    plt.savefig("chi2-test-"+str(lum_arg)+"fb-1_values.pdf")
    plt.cla()
    plt.clf()

    #---- Zscore hist alone with discrete colors
    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
            ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    bounds = np.linspace(0, 5, 6)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Xs), np.log(hist_Q2s), weights=hist_zscores, bins=[
                np.log(Bins['x']), np.log(Bins['Q2'])], cmap=cmap)

    ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])

    cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='uniform', ticks=bounds, boundaries=bounds, format='%1i')

    cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', r'$>$ 5', ' '])

    #cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"${\rm z-score}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$")
    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if hist_CrossSection['non_empty']['total'][iQ2*Bins['Nx']+iX]:
                ax.text(np.log(hist_Xs[iQ2*Bins['Nx']+iX]), np.log(hist_Q2s[iQ2*Bins['Nx']+iX]), '{:.2f}'.format(np.sqrt(hist_chi2s[iQ2*Bins['Nx']+iX])),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=6)

    plt.savefig("Zscore-test-"+str(lum_arg)+"fb-1_values.pdf")
