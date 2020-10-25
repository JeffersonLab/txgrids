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
rc('font', **{'family': 'sans-serif', 'sans-serif': []})
rc('text', usetex=True)
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)

lum_target=100
scenario_choice = "optimistic"
def load_obj(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

def GetSysUnc(sys_path):
    L=open(sys_path).readlines()
    L=[_.strip() for _  in L]
    L=[_ for _ in L if '-----' not in _  if _!='']
    H=L[0].split()
    I=[i for i in range(len(L)) if L[i].startswith('beam')]
    I.append(len(L))

    DATA={}
    for i in range(len(I)-1):
        msg=L[I[i]]
        tab=[[float(v) for v in _.split()] for _ in L[I[i]+1:I[i+1]]]
        tab=np.transpose(tab)
        data={}
        El,Ep=[float(v) for v in msg.split(':')[1].split('x')]
        data['rs']=np.sqrt(4*El*Ep)
        data['x']=tab[0]
        data['Q2']=tab[1]
        data['stat_u(%)']=tab[4]#/np.sqrt(10)
        data['syst_u(%)']=tab[5]
        data['norm_c(%)']=tab[6]
        DATA[msg]=pd.DataFrame(data).query('Q2>1')
    data=pd.concat([DATA[msg] for msg in DATA.keys()], ignore_index=True)

    return data

def GetBins(sys_path):
    #! Barak's
    Bins = {}

    sys = GetSysUnc(sys_path)
    sys = sys.loc[sys['rs'] == rs]

    Bins['x_cv'] = np.array(sorted(set(np.array(sys['x']))))
    Bins['Q2_cv'] = np.array(sorted(set(np.array(sys['Q2']))))

    Bins['x'] = np.zeros(Bins['x_cv'].shape[0]+1)
    Bins["x"][0] = 1e-4 #Bins['x_cv'][0]-(Bins['x_cv'][1]-Bins['x_cv'][0])/2.
    for i,x_cv in enumerate(Bins['x_cv']):
        Bins['x'][i+1] = x_cv+(x_cv-Bins['x'][i])

    Bins['Q2'] = np.zeros(Bins['Q2_cv'].shape[0]+1)
    Bins["Q2"][0] = 1. #Bins['Q2_cv'][0]-(Bins['Q2_cv'][1]-Bins['Q2_cv'][0])/2.
    for i, Q2_cv in enumerate(Bins['Q2_cv']):
        Bins['Q2'][i+1] = Q2_cv+(Q2_cv-Bins['Q2'][i])

    Bins['NQ2'] = len(Bins['Q2'])-1
    Bins['Nx'] = len(Bins['x'])-1

    """#! xiaoxuan (purity based)
    Bins = {}
    kin0, kin = gen_grid(pwd+"../expdata/")
    Bins['x'] = sorted(set(np.array(kin)[:, 0]))
    Bins['Q2'] = sorted(set(np.array(kin)[:, 1]))

    Bins['x_cv'] = sorted(set(np.array(kin0)[:, 0]))[:-1]
    Bins['Q2_cv'] = sorted(set(np.array(kin0)[:, 1]))[:-1]

    Bins['NQ2'] = len(Bins['Q2'])-1
    Bins['Nx'] = len(Bins['x'])-1
    """

    """ #! user defined binning
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
    if unctype == 'sysunc':
        h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_CrossSection[unctype]['uncorr']
                       [key]*100./hist_CrossSection['xsecs'][key], bins=[np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis')
    else:
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


def plt_DetailedZscore(tabnames, Bins, hist_Analysis, non_empty_bins):
    print("--Detailed Zscore histogram plotted--")
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_aspect("equal")

    #fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+tabnames[0].replace("_", "\_") + r'}$, \\' +
    #         r' $\textrm{max: '+tabnames[1].replace("_", "\_") + r'}$, \\' +
    #             r"$N_{rep} = "+str(hist_Analysis['Nrep'])+"$, "+lum_label+r", \textbf{Pessimistic Scenario}", fontsize=10, y=0.98)

    props = dict(boxstyle='square', facecolor='white',
                 edgecolor='gray', alpha=0.5)

    ax.text(0.1, 0.95, r'\hspace{-15pt} \textbf{'+scenario_choice+r' Scenario} \\'+r'$\textrm{($H_0$): '+tabnames[0].replace("_", "\_") + r'}$, \\' + r' $\textrm{($H_1$): '+tabnames[1].replace(
        "_", "\_") + r'}$, \\  $\sqrt{s}=140.7$ GeV, '+lum_label, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)

    ax.text(0.1, 0.7, r'$Z = \left|\frac{\frac{d\sigma^{(H_0)}}{dxdQ^2}-\frac{d\sigma^{(H_1)}}{dxdQ^2}}{\delta^{sys,stat}_{H_0,H_1}}\right|$', transform=ax.transAxes, fontsize=15, verticalalignment='center', bbox=props)


    cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
            ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    bounds = np.linspace(0, 5, 6)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['zscores'], bins=[
        np.log(Bins['x']), np.log(Bins['Q2'])], cmap=cmap, norm=norm)

    ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])

    cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='uniform', ticks=bounds, boundaries=bounds, format='%1i')

    cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', r'$>$ 5', ' '])

    #cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"${\rm Z}$")
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
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.2f}'.format(np.sqrt(
                    hist_Analysis['chi2s'][iQ2*Bins['Nx']+iX])), color="w", ha="center", va="center", fontweight="bold", fontsize=5)
    py.tight_layout()


def plt_DetailedZscore_witherrors(tabnames, Bins, hist_Analysis, non_empty_bins):
    print("--Detailed Zscore with Errors histogram plotted--")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_aspect("equal")

    #pdfbaseline = "NNPDF31_nnlo_pch_as_0118"
    #title = r"{\hspace{12pt} Ultra-"+scenario_choice+r"} \\ ($10\% \times$ systematics)"
    title = r"\textbf{Optimistic}"
    fig.suptitle(r'\hspace{-15pt}\textbf{'+title.replace("_", "\_") + r'}' , fontsize=20, y=0.95)

    #fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+tabnames[0].replace("_", "\_") + r'}$, \\' +
    #         r' $\textrm{max: '+tabnames[1].replace("_", "\_") + r'}$, \\' +
    #             r"$N_{rep} = "+str(hist_Analysis['Nrep'])+"$, "+lum_label+r", \textbf{Pessimistic Scenario}", fontsize=10, y=0.98)

    props = dict(boxstyle='square', facecolor='white',
                 edgecolor='gray', alpha=0.5)

    #ax.text(0.035, 0.97, r'\hspace{-15pt} \textbf{'+scenario_choice+r'}, ' + r" $\sqrt{s}=140.7$ GeV, "+lum_label+r" \\"
    #    +r' $\textrm{\bf($H_0$): '+tabnames[0].replace("_", "\_") + r'}$, ' 
    #    +r' $\textrm{\bf($H_1$): '+tabnames[1].replace( "_", "\_") + r'}$,' \
    #    +r" $N_{rep} = " +str(hist_Analysis['Nrep'])+r'$ \\',\
    #    #r'$ Z_{k\,\,(\text{PDFrep})} = \left|\frac{d\sigma_k^{(H_0)}-d\sigma_k^{(H_1)}}{\delta^{sys,stat}_{H_0,H_1}}\right|$',\
    #    #+r"\\ systematics = "+"{:.2f}".format(hist_Analysis['sys_target'])+r"$\times$original (corr and uncorr)", \
    #    transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)

    #ax.text(0.035, 0.875, r'$ Z_{k\,\,(\text{PDFrep})} = \left|\frac{d\sigma_0^{(H_0)}-d\sigma_k^{(H_1)}}{\delta^{sys,stat}_{H_0,H_1}}\right|$',
    #        transform=ax.transAxes, fontsize=25, verticalalignment='center', bbox=props)

    #ax.text(0.85, 0.1, r'\hspace{-10pt} $\mu_Z $',#\\ \vspace{10pt} $\pm \sigma_Z$',
    #        transform=ax.transAxes, fontsize=25, verticalalignment='center', bbox=props)


    cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
            ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    bounds = np.linspace(0, 5, 6)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['zscores'], bins=[
        np.log(Bins['x']), np.log(Bins['Q2'])], cmap=cmap, norm=norm)

    ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])

    cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='uniform', ticks=bounds, boundaries=bounds, format='%1i')

    cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', r'$>$ 5', ' '])

    #cbar = plt.colorbar()
    #cbar.ax.set_xlabel(r"${\bf \rm <Z>}$", fontsize=15)
    cbar.ax.set_title(r"{\boldmath $Z$}", fontsize=25)
    ax.set_xlabel(r"{\boldmath $x$}", fontsize=25)
    ax.xaxis.set_label_coords(0.95, -0.015)
    ax.set_ylabel(r"{\boldmath $Q^2$}", fontsize=25, rotation=0)
    ax.yaxis.set_label_coords(-0.045, 0.88)
    ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels([r'', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
    ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
    ax.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    ax.tick_params(which='major', direction="in", length=7)
    ax.tick_params(which='minor', direction="in", length=4)

    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if non_empty_bins[iQ2*Bins['Nx']+iX]:
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.1f}'.format(
                    hist_Analysis['zscores_f'][iQ2*Bins['Nx']+iX]),#+"\n"+r'$\pm$'+'{:.2f}'.format(hist_Analysis['std_zscores'][iQ2*Bins['Nx']+iX]),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=9)
    py.tight_layout()

def plt_DetailedPDFweights_witherrors(tabnames, Bins, hist_Analysis, non_empty_bins):
    print("--Detailed Zscore with Errors histogram plotted--")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_aspect("equal")

    #fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+tabnames[0].replace("_", "\_") + r'}$, \\' +
    #         r' $\textrm{max: '+tabnames[1].replace("_", "\_") + r'}$, \\' +
    #             r"$N_{rep} = "+str(hist_Analysis['Nrep'])+"$, "+lum_label+r", \textbf{Pessimistic Scenario}", fontsize=10, y=0.98)

    props = dict(boxstyle='square', facecolor='white',
                 edgecolor='gray', alpha=0.5)

    ax.text(0.05, 0.95, r'\hspace{-15pt} \textbf{'+scenario_choice+r' Scenario} \\'+r'$\textrm{($H_0$): '+tabnames[0].replace("_", "\_") + r'}$, \\' + r' $\textrm{($H_1$): '+tabnames[1].replace(
        "_", "\_") + r'}$, \\' + r"$N_{rep} = "+str(hist_Analysis['Nrep'])+"$, $\sqrt{s}=140.7$ GeV, "+lum_label, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)

    ax.text(0.05, 0.7, r'$w = \sum_k^{N_{rep}}e^{-\frac{1}{2}\chi^2_k} \pm \sigma^{SD}$', transform=ax.transAxes, fontsize=15, verticalalignment='center', bbox=props)


    #cmap = matplotlib.colors.ListedColormap(#['#3535FD', 
    #        ['#2A00D5', '#63009E', '#A1015D', '#D80027', '#FE0002'])

    # define the bins and normalize
    #bounds = np.linspace(0, 5, 6)
    #norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    h = plt.hist2d(np.log(hist_Analysis['x']), np.log(hist_Analysis['Q2']), weights=hist_Analysis['PDFweights'], bins=[
        np.log(Bins['x']), np.log(Bins['Q2'])], cmap='viridis', norm=matplotlib.colors.LogNorm())

    cbar = plt.colorbar()
    #ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])

    #cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
    #    spacing='uniform', ticks=bounds, boundaries=bounds, format='%1i')

    #cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', r'$>$ 5', ' '])

    #cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$w$")
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
                ax.text(np.log(hist_Analysis['x'][iQ2*Bins['Nx']+iX]), np.log(hist_Analysis['Q2'][iQ2*Bins['Nx']+iX]), '{:.1e}'.format(
                    hist_Analysis['PDFweights'][iQ2*Bins['Nx']+iX])+"\n"+r'$\pm$'+'{:.1e}'.format(hist_Analysis['std_PDFweights'][iQ2*Bins['Nx']+iX]),
                        color="w", ha="center", va="center", fontweight="bold", fontsize=4)
    py.tight_layout()

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
    # pwd+"../expdata/data/xQ2binTable-xiaoxuan-060220+syst.npy"
    sys_path = pwd+"../expdata/src/ep_NC_"+scenario_choice+".dat"

    Sample={}
    hist_CrossSection={}
    CovMat={}
    hist_Analysis={}
    #----

    #--physical params and cuts
    #lum = str(lum_arg)+':fb-1'
    lum_label = r'$\mathcal{L} = '+str(lum_target)+r'\,fb^{-1}$'
    cuts_label = r"$ W^2 > 10\,\,,\,\, Q^2> 1$"

    #to tweak:
    #-----------------------------------------------
    #MCEG integration accuracy
    neval = 10000  # default should be: 10000  # (neval/2) ~ batch_size
    nitn = 1000  # default should be: 1000
    min_Q = 0.1

    #SFs
    min_SF = 'Rs=0.5'
    max_SF = 'Rs=1.0'

    seeds = {'min': 3, 'max': 4}
    rs = 140.71247279470288  # 140.7
    #-----------------------------------------------

    #--- acceptance/purity based binning
    Bins = GetBins(sys_path)

    #--- kinematic cuts
    def veto00(x, y, Q2, W2):
        if W2 < 10:
            return 0
        elif Q2 < 1:
            return 0
        elif y > 0.98 or y < 1e-3:
            return 0
        else:
            return 1
    
    if rep != -1:
        #--Getting events
        #if not os.path.isfile(Sample_path):
        #    print(Sample_path+" doesn't exist, please run chi2-test_generate.py first.")
        #    exit(1)
        #else:
        #    Sample = load_obj(Sample_path)

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

        hist_Analysis['sys_target'] = hist_CrossSection['sys_target']
        #--Plotting-----------------------------------------------
        #---------------------------------------------------------
        print("--Plotting--")
        #--General settings
        gs, fig = plt_SetGridSpec((2, 4), (min_SF, max_SF))

        #---- Chi2 hist + sys unc + stat unc 
        title_pos=(0.5,1.04)

        plt_chi2s(0, title_pos, Bins, hist_Analysis)


        #plt_xsecs([0,1],'max', title_pos, Bins, hist_Analysis, Sample)
        #plt_xsecs([1,1],'min', title_pos, Bins, hist_Analysis, Sample)

        plt_unc([0,2],'sysunc','max', title_pos, Bins, hist_Analysis, hist_CrossSection)
        plt_unc([1,2],'sysunc','min', title_pos, Bins, hist_Analysis, hist_CrossSection)

        plt_unc([0,3],'statunc','max', title_pos, Bins, hist_Analysis, hist_CrossSection)
        plt_unc([1,3],'statunc','min', title_pos, Bins, hist_Analysis, hist_CrossSection)

        plt.savefig(sub_wdir+"chi2-overview-"+str(lum_target)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()


        #---- Chi2 hist alone with values written per bin
        print " "
        plt_Detailedchi2s(Bins, hist_Analysis,hist_CrossSection)

        plt.savefig(sub_wdir+"chi2-detailed-"+str(lum_target)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()

        #---- Zscore hist alone with discrete colors
        print " "
        plt_DetailedZscore((min_SF, max_SF), Bins,
                           hist_Analysis, hist_CrossSection['non_empty']['total'])

        plt.savefig(sub_wdir+"Zscore-detailed-"+str(lum_target)+"fb-1_pdfrep"+str(rep)+".pdf")
        plt.cla()
        plt.clf()

    else:
        ##------------
        Nrep =0
        (_, dirnames, filenames) = walk(sub_wdir).next()
        tot_hist_Analysis = {}
        tot_hist_Analysis['chi2s'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['zscores'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['PDFweights'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['std_chi2s'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['std_zscores'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['std_PDFweights'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['x'] = np.zeros(Bins['Nx']*Bins['NQ2'])
        tot_hist_Analysis['Q2'] = np.zeros(Bins['Nx']*Bins['NQ2'])

        non_empty_bins = np.zeros(Bins['Nx']*Bins['NQ2']) > 1
        ONCE=True
        for dirname in dirnames:
            if os.path.isfile(sub_wdir+dirname+"/hist_Analysis.p"):
                hist_Analysis = load_obj(sub_wdir+dirname+"/hist_Analysis.p")
                if ONCE:
                    hist_CrossSection = load_obj(sub_wdir+dirname+"/hist_CrossSection.p")
                    tot_hist_Analysis['sys_target'] = hist_CrossSection['sys_target']
                    ONCE=False

                tot_hist_Analysis['chi2s'] += hist_Analysis['chi2s']
                tot_hist_Analysis['zscores'] += np.sqrt(hist_Analysis['chi2s'])
                tot_hist_Analysis['PDFweights'] += hist_Analysis['PDFweights']

                hist_CrossSection = load_obj(sub_wdir+dirname+"/hist_CrossSection.p")
                non_empty_bins += hist_CrossSection['non_empty']['total']
                Nrep+=1

        tot_hist_Analysis['chi2s'] /= Nrep
        tot_hist_Analysis['zscores'] /= Nrep
        tot_hist_Analysis['PDFweights'] /= Nrep
        tot_hist_Analysis['zscores_f'] = np.copy(tot_hist_Analysis['zscores'])
        tot_hist_Analysis['Nrep'] = Nrep

        for dirname in dirnames:
            if os.path.isfile(sub_wdir+dirname+"/hist_Analysis.p"):
                hist_Analysis = load_obj(sub_wdir+dirname+"/hist_Analysis.p")
                tot_hist_Analysis['std_chi2s'] += (hist_Analysis['chi2s']-tot_hist_Analysis['chi2s'])**2
                tot_hist_Analysis['std_zscores'] += (np.sqrt(hist_Analysis['chi2s'])-tot_hist_Analysis['zscores'])**2
                tot_hist_Analysis['std_PDFweights'] += (hist_Analysis['PDFweights']-tot_hist_Analysis['PDFweights'])**2

        tot_hist_Analysis['std_chi2s'] /= Nrep
        tot_hist_Analysis['std_chi2s'] = np.sqrt(tot_hist_Analysis['std_chi2s'])

        tot_hist_Analysis['std_zscores'] /= Nrep
        tot_hist_Analysis['std_zscores'] = np.sqrt(tot_hist_Analysis['std_zscores'])

        tot_hist_Analysis['std_PDFweights'] /= Nrep
        tot_hist_Analysis['std_PDFweights'] = np.sqrt(tot_hist_Analysis['std_PDFweights'])


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
        tot_hist_Analysis['PDFweights'][tot_hist_Analysis['PDFweights'] == 0] = np.nan

        tot_hist_Analysis['zscores'][tot_hist_Analysis['zscores'] == -1] = 0

        plt_DetailedZscore_witherrors((min_SF,max_SF), Bins, tot_hist_Analysis, non_empty_bins)

        plt.savefig(sub_wdir+"ZscoreWithError-detailed-" +
                    str(lum_target)+"fb-1_Nrep"+str(Nrep)+".pdf")
        plt.cla()
        plt.clf()

        plt_DetailedPDFweights_witherrors((min_SF,max_SF), Bins, tot_hist_Analysis, non_empty_bins)

        plt.savefig(sub_wdir+"PDFweightsWithError-detailed-" +
                    str(lum_target)+"fb-1_Nrep"+str(Nrep)+".pdf")
        plt.cla()
        plt.clf()
