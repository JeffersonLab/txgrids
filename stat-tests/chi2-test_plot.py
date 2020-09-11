'''
 # @ Author: Rabah Abdul Khalek
 # @ Create Time: 2020-04-24 14:02:30
 # @ Modified by: Rabah Abdul Khalek
 # @ Modified time: 2020-04-24 15:04:17
 # @ Description: New ideas for impact studies at  EIC
 '''

import sys,os
from os import walk
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory.tools import save, load, convert_lum
from theory.mceg  import MCEG
from expdata.driver import gen_grid
import pickle

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

def load_obj(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

def GetBins(pwd=""):
    Bins = {}
    kin0, kin = gen_grid(pwd+"../expdata/")
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

def units(units):
    one=0.3893793721 #--GeV2 mbarn
    if   units=='GeV^-2':return 1
    elif units=='fb'    :return one*1e12 
    else: sys.exit('%s conversion not available')

def plt_SetGridSpec(size,tabnames):
    nrows, ncols = size
    widths = [2, 1, 1, 1]
    heights = [1, 1]
    gs = gridspec.GridSpec(nrows, ncols, width_ratios=widths,
                           height_ratios=heights)
    fig = py.figure(figsize=(ncols*5, nrows*3.5))  # , constrained_layout=True)

    fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+tabnames[0].replace("_", "\_") + r'}$,' +
                 r' $\textrm{max: '+tabnames[1].replace("_", "\_") + r'}$, ' +
                 lum_label, fontsize=10, y=0.98)

    fig.subplots_adjust(top=0.85, bottom=0.15)

    return gs,fig

def plt_chi2s(gs_ind, title_pos, Bins, hist_Analysis):
    print("--chi2 histogram plotted--")

    title_xpos, title_ypos = title_pos

    #---- chi2 hist
    ax = py.subplot(gs[:, gs_ind])
    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['chi2s'], bins=[
        np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())

    ax.text(title_xpos, title_ypos+0.02, r'$\chi^2_{tot}/N_{bins}$ = '+' %0.4e' % hist_Analysis['tot_chi2'],
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=15)

    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\chi^2_{bin}$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$',
                        r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

def plt_xsecs(gs_ind, key, title_pos, Bins, hist_Analysis, Sample):

    print("--xsec["+key+"] histogram plotted--")

    title_xpos, title_ypos = title_pos
    ax = py.subplot(gs[gs_ind[0], gs_ind[1]])
    ax.text(title_xpos, title_ypos, r'$\sigma^{'+key+'}_{tot}$'+' = %0.2e' % (Sample[key]['tot_xsec']) + r'$\pm$'+'%0.2e' % (Sample[key]['var_xsec'])+r' [fb]',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_CrossSection['xsecs'][key], bins=[
        np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())

    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\frac{d\sigma_{"+key+"}}{dxdQ^2}$")

    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    #ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$',
                        r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

def plt_unc(gs_ind, unctype, key, title_pos, Bins, hist_Analysis, hist_CrossSection):

    print("--"+unctype+"["+key+"] histogram plotted--")
    title_xpos, title_ypos = title_pos

    uncname = unctype.replace('unc', '')

    ax = py.subplot(gs[gs_ind[0], gs_ind[1]])
    ax.text(title_xpos, title_ypos, r''+uncname+' Unc ('+key+')',  # r'$\delta^{max,sys}_{tot}$'+' = %0.2e' %
            #((np.sqrt(Sample['max']['Nevents'])/lum_arg)
            # * 100./Sample['max']['tot_xsec']) + r' \%',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10)
    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_CrossSection[unctype]
                   [key]*100./hist_CrossSection['xsecs'][key], bins=[np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\delta^{"+key+"}_{"+uncname+"} [\%]$")
    #cbar.ax.set_xlabel(r"$\frac{1}{\sigma_{tot}}\frac{d\sigma_{max}}{dxdQ^2}$")
    #ax.set_xlabel(r"$x$")

    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-4}$', r'$10^{-3}$',
                        r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)


def plt_Detailedchi2s(Bins, hist_Analysis, hist_CrossSection):
    print("--Detailed chi2 histogram plotted--")
    fig, ax = plt.subplots()
    ax.set_aspect("equal")


    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['chi2s'], bins=[
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
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.2f}'.format(hist_Analysis['chi2s'][iQ2*Bins['Nx']+iX]),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=6)


def plt_DetailedZscore(Bins, hist_Analysis, hist_CrossSection):
    print("--Detailed Zscore histogram plotted--")
    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
            ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    bounds = np.linspace(0, 5, 6)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['zscores'], bins=[
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
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.2f}'.format(np.sqrt(hist_Analysis['chi2s'][iQ2*Bins['Nx']+iX])),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=6)


def plt_DetailedZscore_witherrors(Bins, hist_Analysis, non_empty_bins):
    print("--Detailed Zscore with Errors histogram plotted--")
    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
            ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    bounds = np.linspace(0, 5, 6)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['zscores'], bins=[
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
            if non_empty_bins[iQ2*Bins['Nx']+iX]:
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.2f}'.format(
                    hist_Analysis['zscores_f'][iQ2*Bins['Nx']+iX])+"\n"+r'$\pm$'+'{:.1f}'.format(hist_Analysis['std_zscores'][iQ2*Bins['Nx']+iX]),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=5)

if __name__ == "__main__":

    #--- Getting arguments from user
    if len(sys.argv) < 2:
        print "usage: ./chi2-test_plot.py <Luminosity (fb^-1)> <pdf_replica=-1 for total analysis>"
        exit(1)

    lum_arg = float(sys.argv[1])
    rep = int(sys.argv[2])
    #---

    #--- handling paths
    pwd = os.getcwd()+"/"

    wdir = pwd+"output/"
    if not os.path.isdir(wdir): os.system("mkdir "+wdir)

    sub_wdir = wdir+"lum"+str(lum_arg)+"/"
    if not os.path.isdir(sub_wdir): os.system("mkdir "+sub_wdir)

    sub_sub_wdir = ""
    if rep != -1:
        sub_sub_wdir = sub_wdir+"pdfrep" + str(rep)+"/"
        if not os.path.isdir(sub_sub_wdir): os.system("mkdir "+sub_sub_wdir)

    Sample_path = sub_sub_wdir+"Sample.p"
    hist_CrossSection_path = sub_sub_wdir+"hist_CrossSection.p"
    CovMat_path = sub_sub_wdir+"CovMat.p"
    hist_Analysis_path = sub_sub_wdir+"hist_Analysis.p"

    Sample={}
    hist_CrossSection={}
    CovMat={}
    hist_Analysis={}
    #----

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
    Bins=GetBins(pwd)

    #--- kinematic cuts
    def veto00(x,y,Q2,W2):
        if   W2 < 10       : return 0
        elif Q2 < 1        : return 0
        else               : return 1
    
    if rep != -1:
        #--Getting events
        if not os.path.isfile(Sample_path):
            print(Sample_path+" doesn't exist, please run chi2-test_generate.py first.")
            exit(1)
        else:
            Sample = load_obj(Sample_path)

        #--Getting CrossSections
        if not os.path.isfile(hist_CrossSection_path):
            print(hist_CrossSection_path+" doesn't exist, please run chi2-test_generate.py first.")
            exit(1)
        else:
            hist_CrossSection = load_obj(hist_CrossSection_path)

        #--Getting Analysis: Chi2s and Z-scores
        if not os.path.isfile(hist_Analysis_path):
            print(hist_Analysis_path+" doesn't exist, please run chi2-test_generate.py first.")
            exit(1)
        else:
            hist_Analysis = load_obj(hist_Analysis_path)

        #--Plotting-----------------------------------------------
        #---------------------------------------------------------
        print("--Plotting--")

        #--General settings
        gs, fig = plt_SetGridSpec((2, 4), (min_SF, max_SF))

        #---- Chi2 hist + sys unc + stat unc 
        title_pos=(0.5,1.04)

        plt_chi2s(0, title_pos, Bins, hist_Analysis)


        plt_xsecs([0,1],'max', title_pos, Bins, hist_Analysis, Sample)
        plt_xsecs([1,1],'min', title_pos, Bins, hist_Analysis, Sample)

        plt_unc([0,2],'sysunc','max', title_pos, Bins, hist_Analysis, hist_CrossSection)
        plt_unc([1,2],'sysunc','min', title_pos, Bins, hist_Analysis, hist_CrossSection)

        plt_unc([0,3],'statunc','max', title_pos, Bins, hist_Analysis, hist_CrossSection)
        plt_unc([1,3],'statunc','min', title_pos, Bins, hist_Analysis, hist_CrossSection)

        plt.savefig(sub_wdir+"chi2-overview-"+str(lum_arg)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()


        #---- Chi2 hist alone with values written per bin
        print " "
        plt_Detailedchi2s(Bins, hist_Analysis,hist_CrossSection)

        plt.savefig(sub_wdir+"chi2-detailed-"+str(lum_arg)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()

        #---- Zscore hist alone with discrete colors
        print " "
        plt_DetailedZscore(Bins,hist_Analysis,hist_CrossSection)

        plt.savefig(sub_wdir+"Zscore-detailed-"+str(lum_arg)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()

    else:
        (_, dirnames, filenames) = walk(sub_wdir).next()
        tot_hist_Analysis = {}
        tot_hist_Analysis['chi2s'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['zscores'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['std_chi2s'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['std_zscores'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['x'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['Q2'] = np.zeros(Bins['Nx']*Bins['NQ2'])

        non_empty_bins = np.zeros(Bins['Nx']*Bins['NQ2']) > 1

        for dirname in dirnames:
            hist_Analysis = load_obj(sub_wdir+"/"+dirname+"/hist_Analysis.p")
            tot_hist_Analysis['chi2s'] += hist_Analysis['chi2s']
            tot_hist_Analysis['zscores'] += np.sqrt(hist_Analysis['chi2s'])

            hist_CrossSection = load_obj(sub_wdir+"/"+dirname+"/hist_CrossSection.p")
            non_empty_bins += hist_CrossSection['non_empty']['total']

        tot_hist_Analysis['chi2s'] /= len(dirnames)
        tot_hist_Analysis['zscores'] /= len(dirnames)
        tot_hist_Analysis['zscores_f'] = np.copy(tot_hist_Analysis['zscores'])

        for dirname in dirnames:
            hist_Analysis = load_obj(sub_wdir+"/"+dirname+"/hist_Analysis.p")
            tot_hist_Analysis['std_chi2s'] += (hist_Analysis['chi2s']-tot_hist_Analysis['chi2s'])**2
            tot_hist_Analysis['std_zscores'] += (np.sqrt(hist_Analysis['chi2s'])-tot_hist_Analysis['zscores'])**2

        tot_hist_Analysis['std_chi2s'] /= len(dirnames)
        tot_hist_Analysis['std_chi2s'] = np.sqrt(tot_hist_Analysis['std_chi2s'])

        tot_hist_Analysis['std_zscores'] /= len(dirnames)
        tot_hist_Analysis['std_zscores'] = np.sqrt(tot_hist_Analysis['std_zscores'])

        #--- for plotting purposes
        for iQ2 in range(0, Bins['NQ2']):
            for iX in range(0, Bins['Nx']):
                if tot_hist_Analysis['zscores'][iQ2*Bins['Nx']+iX]<1:
                    tot_hist_Analysis['zscores'][iQ2*Bins['Nx']+iX]=-1
                else:
                    if not np.isnan(tot_hist_Analysis['zscores'][iQ2*Bins['Nx']+iX]):
                        tot_hist_Analysis['zscores'][iQ2*Bins['Nx']+iX] = int(tot_hist_Analysis['zscores'][iQ2*Bins['Nx']+iX])

                tot_hist_Analysis['x'][iQ2*Bins['Nx']+iX] = Bins['x_cv'][iX]
                tot_hist_Analysis['Q2'][iQ2*Bins['Nx']+iX] = Bins['Q2_cv'][iQ2]

        tot_hist_Analysis['zscores'][tot_hist_Analysis['zscores'] >= 5] = 5

        tot_hist_Analysis['chi2s'][tot_hist_Analysis['chi2s']==0] = np.nan
        tot_hist_Analysis['zscores'][tot_hist_Analysis['zscores'] == 0] = np.nan

        tot_hist_Analysis['zscores'][tot_hist_Analysis['zscores'] == -1] = 0

        plt_DetailedZscore_witherrors(Bins, tot_hist_Analysis, non_empty_bins)

        plt.savefig(sub_wdir+"ZscoreWithError-detailed-" +
                    str(lum_arg)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()
