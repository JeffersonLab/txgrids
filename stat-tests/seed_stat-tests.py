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
from theory.tools import save, load
from theory.seed_mceg  import MCEG
from ndtest import ndtest

#--matplotlib
import matplotlib
import matplotlib.pyplot as plt
import pylab  as py
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

from scipy import stats
import pandas as pd
from collections import OrderedDict, Counter

from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

wdir='.stats-tests'


#to tweak:
#-----------------------------------------------
if len(sys.argv) < 1:
    print "usage: ./stat-tests.py <Nevents>"
    exit(1)

Nevents = int(sys.argv[1])

#--physical params and cuts
rs= 140.7
lum='10:fb-1'
lum_label= r'$\mathcal{L} = 10\,fb^{-1}$'
sign=1 #--electron=1 positron=-1

cuts_label = r"$ W^2 > 10\,\,,\,\, Q^2> 1$"
def veto00(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1


# "Kolmogorov-Smirnov"]  # "Mann-Whitney",
stat_methods = ["Kolmogorov-Smirnov"]
                #"Kruskal-Wallis"]  # "Barlett" #"Fligner-Killeen" #"Ansari-Bradley"
fig_composition = [3,2] # make room for 2 additional hists
#-----------------------------------------------

#--lhapdf set and stf idx
s = {'min':{}, 'max':{}}
"""
s['min']['tabname'] = 'JAM4EIC'
s['min']['iset'] = 0
s['min']['iF2'] = 90001
s['min']['iFL'] = 90002
s['min']['iF3'] = 90003

"""
s['min']['tabname'] = 'NNPDF31_nnlo_pch_as_0118_SF'
s['min']['iset'] = 0
s['min']['iF2'] = 908
s['min']['iFL'] = 909
s['min']['iF3'] = 910
s['min']['seed'] = 1


s['max']['tabname'] = 'NNPDF31_nnlo_pch_as_0118_SF'
s['max']['iset'] = 0
s['max']['iF2'] = 908
s['max']['iFL'] = 909
s['max']['iF3'] = 910
s['max']['seed'] = 1234567

#common keys
for key in s.keys():
    s[key]['wdir']      = wdir
    s[key]['sign']      = sign
    s[key]['rs'] = rs
    s[key]['fname'] = 'mceg00'
    s[key]['veto'] = veto00
    s[key]['fdata'] = wdir+"/"+s[key]['tabname']+"_data"+str(int(Nevents/1000))+"k_"+str(s[key]['seed'])+".po"


#--generate events
#output keys: ['Y', 'X', 'Q2', 'W', 'rs']
for key in s.keys():
    if not os.path.isfile(s[key]['fdata']):
        mceg=MCEG(**(s[key]))
        mceg.buil_mceg()
        s[key].update(mceg.gen_events(Nevents))
        save(s[key], s[key]['fdata'])
    else:
        s[key] = load(s[key]['fdata'])

np.set_printoptions(threshold=sys.maxsize)

smin=np.array([s['min']['W'],s['min']['X'],s['min']['Q2']])
smax=np.array([s['max']['W'],s['max']['X'],s['max']['Q2']])

"""
#------ 
#--------------------------------------------------------------------------------------------------------------------
stat, pvalue = stats.ks_2samp(smin[0,:], smax[0,:])
print "Kolmogorov-Smirnov: stat = ", stat, " p-value = ", pvalue, "..."
print "----------"
#--------------------------------------------------------------------------------------------------------------------
"""
#------
#--------------------------------------------------------------------------------------------------------------------

#to tweak:
#---
NQ2bins = 10 # int(np.log10(Nevents)*25)
print "NQ2bins = ", NQ2bins
NXbins = 10 # int(np.log10(Nevents)*25)
print "NXbins = ", NXbins
choice_bins = 'min'
#---

Q2bins_smin = np.logspace(np.log10(np.min(s['min']['Q2'])), np.log10(np.max(s['min']['Q2'])),NQ2bins+1)
Q2bins_smax = np.logspace(np.log10(np.min(s['max']['Q2'])), np.log10(np.max(s['max']['Q2'])),NQ2bins+1)

Xbins_smin = np.logspace(np.log10(np.min(s['min']['X'])), np.log10(np.max(s['min']['X'])),NXbins+1)
Xbins_smax = np.logspace(np.log10(np.min(s['max']['X'])), np.log10(np.max(s['max']['X'])),NXbins+1)

if choice_bins == 'min':
    Q2bins = Q2bins_smin
    Xbins = Xbins_smin
    resultpath = 'plots/ks-test_perxQ2_'+str(NXbins)+"-"+str(NQ2bins)+'sminbins_'+str(Nevents/1000)+'k-events.pdf'
elif choice_bins == 'max':
    Q2bins = Q2bins_smax
    Xbins = Xbins_smax
    resultpath = 'plots/ks-test_perxQ2_'+str(NXbins)+"-"+str(NQ2bins)+'smaxbins_'+str(Nevents/1000)+'k-events.pdf'

smin_perQ2bin_perXbin = []
smax_perQ2bin_perXbin = []

Q2_dig_smin = np.digitize(s['min']['Q2'], Q2bins)
Q2_dig_smax = np.digitize(s['max']['Q2'], Q2bins)

X_dig_smin = np.digitize(s['min']['X'], Xbins)
X_dig_smax = np.digitize(s['max']['X'], Xbins)

empty_bins=[]
for iQ2 in range(1, NQ2bins+1):
    Q2mask_smin = np.where(Q2_dig_smin == iQ2, True, False)
    Q2mask_smax = np.where(Q2_dig_smax == iQ2, True, False)
    smin_perQ2bin_perXbin.append([])
    smax_perQ2bin_perXbin.append([])
    for iX in range(1, NXbins+1):
        Xmask_smin = np.where(X_dig_smin == iX, True, False)
        Xmask_smax = np.where(X_dig_smax == iX, True, False)

        #combining masks
        mask_smin = Q2mask_smin * Xmask_smin
        mask_smax = Q2mask_smax * Xmask_smax

        #apply mask
        masked_smin = smin[:, mask_smin]
        masked_smax = smax[:, mask_smax]

        #append
        smin_perQ2bin_perXbin[iQ2-1].append(masked_smin)
        smax_perQ2bin_perXbin[iQ2-1].append(masked_smax)

        #keep track of empty bins
        if not list(masked_smin[0]) or not list(masked_smax[0]): #! TBD: if this should be an 'and' instead
            empty_bins.append([iQ2-1,iX-1])

fig = plt.figure()
fig.subplots_adjust(top=0.85, bottom=0.)
fig.suptitle(r'\hspace{-15pt}$\textrm{min: '+s['min']['tabname'].replace("_", "\_") + r'}$' +
             r'\\ $\textrm{max: '+s['max']['tabname'].replace("_", "\_") + r'}$', fontsize=10, y=0.98)
for imethod, stat_method in enumerate(stat_methods): # ["Kolmogorov-Smirnov", "Ansari-Bradley", "Barlett", "Fligner-Killeen", "Kruskal-Wallis"]
    print "method: ", stat_method, "... "
    pvalues_perQ2bin_perXbin = []
    for iQ2 in range(0, NQ2bins):
        pvalues_perQ2bin_perXbin.append([])
        for iX in range(0, NXbins):
            if [iQ2, iX] in empty_bins:
                pvalue=1
            else:
                if stat_method == "Kolmogorov-Smirnov":
                    stat, pvalue = stats.ks_2samp(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
                elif stat_method == "Ansari-Bradley":
                    stat, pvalue = stats.ansari(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
                elif stat_method == "Mann-Whitney":
                    stat, pvalue = stats.mannwhitneyu(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
                elif stat_method == "Barlett":
                    stat, pvalue = stats.bartlett(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
                elif stat_method == "Fligner-Killeen":
                    stat, pvalue = stats.fligner(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
                elif stat_method == "Kruskal-Wallis":
                    stat, pvalue = stats.kruskal(smin_perQ2bin_perXbin[iQ2][iX][0], smax_perQ2bin_perXbin[iQ2][iX][0])
            
            pvalues_perQ2bin_perXbin[iQ2].append(pvalue)

    #plt.clf()
    #fig = plt.figure()
    #fig.suptitle(r'$\textrm{'+stat_method+r'\,\,test}$' +
    #            r'\\ $\textrm{min: '+s['min']['tabname'].replace("_", "\_") + r'}$'+
    #            r'\\ $\textrm{max: '+s['max']['tabname'].replace("_", "\_") + r'}$', fontsize=10)

    ax = fig.add_subplot(fig_composition[0]*100+fig_composition[1]*10+imethod+1)
    plt.tight_layout(rect=[0, 0.08, 1., 0.95])
    ax.title.set_text(r'$\textrm{'+stat_method+r'\,\,test}$')# +
                      #r'\\ $\textrm{min: '+s['min']['tabname'].replace("_", "\_") + r'}$' +
                      #r'\\ $\textrm{max: '+s['max']['tabname'].replace("_", "\_") + r'}$')#, fontsize=10)

    central_Q2bins = np.array(Q2bins)[:-1]+(np.array(Q2bins)[1:]-np.array(Q2bins)[:-1])/2
    central_Xbins = np.array(Xbins)[:-1]+(np.array(Xbins)[1:]-np.array(Xbins)[:-1])/2

    hist_Q2bins = []
    hist_Xbins = []
    hist_pvalues = []
    hist_Nevents = {'max':[],'min':[]}
    Nbins=0
    #reshaping
    for iQ2 in range(0, NQ2bins):
        for iX in range(0, NXbins):
            #if [iQ2, iX] not in empty_bins:
            hist_Xbins.append(central_Xbins[iX])
            hist_Q2bins.append(central_Q2bins[iQ2])
            hist_pvalues.append(pvalues_perQ2bin_perXbin[iQ2][iX])
            hist_Nevents['min'].append(len(list(smin_perQ2bin_perXbin[iQ2][iX][0])))
            hist_Nevents['max'].append(len(list(smax_perQ2bin_perXbin[iQ2][iX][0])))

    #convert p-values to Z-scores
    hist_Z_score=[]
    for i in range(0, NQ2bins*NXbins):
            if hist_pvalues[i] <= 1e-15:
                sigma_lvl = stats.norm.ppf(1-(1e-15)/2)
            elif hist_pvalues[i] == 1:
                sigma_lvl=0
            else:
                sigma_lvl = stats.norm.ppf(1-(hist_pvalues[i])/2)

            hist_Z_score.append(sigma_lvl)

    #hist_Z_score = stats.norm.ppf(1-(np.array(hist_pvalues))/2)

    h = plt.hist2d(np.log(hist_Xbins), np.log(hist_Q2bins), weights=hist_Z_score, bins=[
                NXbins, NQ2bins], cmap='hot_r')
    #plt.clim(0, 1)
    plt.clim(0,5)
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$\sigma-level$")  # , rotation=270)
    #cbar.ax.set_yscale('linear')  # , rotation=270)

    if (imethod+1)%2 == 0:
        #remove xticks labels
        ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
        ax.set_xticklabels([])
        ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
        ax.set_yticklabels([])
    else:
        #remove xticks labels
        ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
        ax.set_xticklabels([])
        ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
        ax.set_yticklabels([r'$1$', r'$10$', r'$100$', r'$1000$', r'$10000$'])
        ax.set_ylabel(r'$Q^2$', size=20)

#ax.text(0.1, 0.7, r'$N_x = %s$, $N_{Q^2} = %s$ \\ $N_{events}= %s$ \\$\sqrt{s}=%0.2f{\rm~GeV}$' %(NXbins, NQ2bins, Nevents, rs), transform=ax.transAxes, size=15)
props = dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.5)
ax.text(1.45, 0.9, r'\hspace{-15pt}$N_x = %s$, $N_{Q^2} = %s$ \\ $N_{events}= %sk$ \\$\sqrt{s}=%0.2f{\rm~GeV}$ \\ %s' % (
    NXbins, NQ2bins, Nevents/1000, rs, cuts_label), transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)

#plot number of events in bins
for ikey,key in enumerate(s.keys()):
    print "Nevents of sample ",key," plotted..."
    ax = fig.add_subplot(fig_composition[0]*100+fig_composition[1]*10+5+ikey)
    ax.title.set_text(r'$\textrm{'+key+r' sample}\,\,N_{events}$')  # +
    h = plt.hist2d(np.log(hist_Xbins), np.log(hist_Q2bins), weights=hist_Nevents[key], bins=[
        NXbins, NQ2bins], cmap='hot_r')
    
    #plt.clim(0, 5)
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(r"$N_{events}$")  # , rotation=270)

    if (ikey+1) % 2 == 0:
        ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
        ax.set_xticklabels([r'$0.0001$', r'$0.001$', r'$0.01$', r'$0.1$'])
        ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
        ax.set_yticklabels([])
        ax.set_xlabel(r'$x$', size=20)

    else:
        ax.set_xticks(np.log([1e-4, 1e-3, 1e-2, 1e-1]))
        ax.set_xticklabels([r'$0.0001$', r'$0.001$', r'$0.01$', r'$0.1$'])
        ax.set_yticks(np.log([1, 10, 100, 1000, 10000]))
        ax.set_yticklabels([r'$1$', r'$10$', r'$100$', r'$1000$', r'$10000$'])
        ax.set_ylabel(r'$Q^2$', size=20)
        ax.set_xlabel(r'$x$', size=20)


"""
plt.hist(central_bins, weights=pvalues_perQ2bin, bins=Q2bins,
        histtype="step", color="blue", ls="solid", linewidth=1.0, label=r"p-value")
plt.xscale('log')
plt.gcf().subplots_adjust(bottom=0.15)
plt.xlabel(r"$Q^2$",fontsize=13)
plt.ylabel(r"$p-value$",fontsize=17)
plt.ylim(0,1.)
plt.title("min: "+s['min']['tabname'].replace("_", "\_") +
        "\n max: "+s['max']['tabname'].replace("_", "\_"), fontsize=12)
handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor(), ls=h.get_linestyle()) for h in handles]
#zer = np.zeros(NQ2bins)
#ax.plot(central_bins, zer, linewidth=1, ls='--', color='black')
#plt.legend(handles=new_handles, labels=labels, fontsize=13, loc='best')

props = dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.5)
ax.text(0.7, 0.95, r"$N^{Q^2}_{bins}=$"+str(NQ2bins)+"\n"+"Nevents="+str(Nevents/1000)+"k", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
"""
plt.savefig(resultpath)
plt.clf
print(resultpath+" saved...")
print "----------"
#--------------------------------------------------------------------------------------------------------------------


"""

#------ Pearson correlation coefficient
#--------------------------------------------------------------------------------------------------------------------
#       Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html#scipy.stats.pearsonr
#!      Notes: - The calculation of the p-value relies on the assumption that each dataset is normally distributed
#P_corr = stats.pearsonr(smin[0],smax[0])
#print P_corr
#--------------------------------------------------------------------------------------------------------------------



#------ Spearman correlation coefficient
#--------------------------------------------------------------------------------------------------------------------
#       Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html#scipy.stats.spearmanr
#!      Notes:  - does not assume that both datasets are normally distributed.
#!              - measure of the monotonicity of the relationship between two datasets.
#example1
#np.random.seed(1234321)
#x2n = np.random.randn(2, 100)
#y2n = np.random.randn(2, 100)
#S_corr = stats.spearmanr(x2n, y2n, axis=1)

#example2
#x1 = np.linspace(1, 100, 100)
#x2 = np.linspace(100,1,100)
#x1 = np.array([x1,x2])
#rho, pvalue = stats.spearmanr(x1, x1, axis=1)
#print rho
#print pvalue

rho, pvalue = stats.spearmanr(smin,smax,axis=1)
columns = ['W_min', 'X_min', 'Q2_min', 'W_max', 'X_max', 'Q2_max']

rho_output = OrderedDict()
for i,key in enumerate(columns):
    rho_output[key]=rho[i]
rho_output = pd.DataFrame.from_dict( rho_output, orient='index', columns=columns)
print(rho_output)

pvalue_output = OrderedDict()
for i,key in enumerate(columns):
    pvalue_output[key]=pvalue[i]
pvalue_output = pd.DataFrame.from_dict( pvalue_output, orient='index', columns=columns)
print(pvalue_output)
print "----------"
#--------------------------------------------------------------------------------------------------------------------



#------
#--------------------------------------------------------------------------------------------------------------------
stat, pvalue = stats.ansari(smin[0, :], smax[0, :])
print "Ansari-Bradley: stat = ", stat, " p-value = ", pvalue, "..."
print "----------"
#--------------------------------------------------------------------------------------------------------------------

#------
#--------------------------------------------------------------------------------------------------------------------
stat, pvalue = stats.bartlett(smin[0, :], smax[0, :])
print "Barlett: stat = ", stat, " p-value = ", pvalue, "..."
print "----------"
#--------------------------------------------------------------------------------------------------------------------

#------
#--------------------------------------------------------------------------------------------------------------------
stat, critical_values, significance_level = stats.anderson_ksamp([smin[0, :], smax[0, :]])
print "Anderson-Darling: stat = ", stat, " critical values = ", pvalue, " significance level = ", significance_level, "..."
print "----------"
#--------------------------------------------------------------------------------------------------------------------

#------
#--------------------------------------------------------------------------------------------------------------------
stat, pvalue = stats.fligner(smin[0, :], smax[0, :])
print "Fligner-Killeen: stat = ", stat, " p-value = ", pvalue, "..."
print "----------"
#--------------------------------------------------------------------------------------------------------------------


#------
#--------------------------------------------------------------------------------------------------------------------
stat, pvalue = stats.kruskal(smin[0:], smax[0:])
print "Kruskal-Wallis: stat = ", stat, " p-value = ", pvalue, "..." 
print "----------"
#--------------------------------------------------------------------------------------------------------------------
"""
os.system("rm -r bins")
os.mkdir('bins/')
i=0
for iQ2 in range(0, NQ2bins):
    for iX in range(0, NXbins):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #plt.tight_layout(rect=[0, 0.08, 1., 0.95])
        ax.title.set_text(r'$\sigma-lvl=$' +
                          str(round(hist_Z_score[i], 2))+"; x="+str(round(central_Xbins[iX], 2))+"; Q2="+str(round(central_Q2bins[iQ2], 2)))

        (n, bins, patches) = ax.hist(smin_perQ2bin_perXbin[iQ2][iX][0], bins='fd', facecolor='red', histtype="step", alpha=1., label='min', ls='-', linewidth=1.5)
        ax.hist(smax_perQ2bin_perXbin[iQ2][iX][0], bins=bins, facecolor='blue',
                histtype="step", alpha=1., label='min', ls='-', linewidth=1.5)

        ax.set_xscale('log')
        ax.set_yscale('log')

        plt.savefig("bins/Q2_"+str(iQ2)+"_X_"+str(iX)+"_p_"+str(round(hist_Z_score[i],2))+".pdf")
        plt.clf
        i+=1

