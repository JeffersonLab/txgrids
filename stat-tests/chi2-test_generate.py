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

def save_obj(obj, path ):
    with open(path, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

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

def GetEvents(path,lum_arg,veto,sign,rep=0):
    """
    sign: --electron=1 positron=-1

    """
    Sample = {'min':{}, 'max':{}}

    Sample['min']['tabname'] = min_SF
    Sample['min']['iset'] = rep
    Sample['min']['iF2'] = 908
    Sample['min']['iFL'] = 909
    Sample['min']['iF3'] = 910

    Sample['max']['tabname'] = max_SF
    Sample['max']['iset'] = rep
    Sample['max']['iF2'] = 908
    Sample['max']['iFL'] = 909
    Sample['max']['iF3'] = 910

    #common keys
    for key in Sample.keys():
        Sample[key]['wdir'] = path
        Sample[key]['sign']      = sign
        Sample[key]['rs'] = rs
        Sample[key]['fname'] = 'mceg00'
        Sample[key]['veto'] = veto
        Sample[key]['fdata'] = path

    #--generate events
    #output keys: ['Y', 'X', 'Q2', 'W', 'rs']
    print("--Getting Events--")
    for key in Sample.keys():
        np.random.seed(seeds[key])
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

    keys_of_interest = ['W','X','Q2']

    for key1 in Sample.keys():
        for key2 in keys_of_interest:
            Sample[key1][key2]=np.array(Sample[key1][key2])

    return Sample

def GetBins(sys_path):
    #! Barak's
    Bins = {}

    sys = GetSysUnc(sys_path)
    sys = sys.loc[sys['rs'] == rs]

    Bins['x_cv'] = np.array(sorted(set(np.array(sys['x']))))
    Bins['Q2_cv'] = np.array(sorted(set(np.array(sys['Q2']))))

    Bins['x'] = np.zeros(Bins['x_cv'].shape[0]+1)
    Bins["x"][0] = Bins['x_cv'][0]-(Bins['x_cv'][1]-Bins['x_cv'][0])/2.
    for i,x_cv in enumerate(Bins['x_cv']):
        Bins['x'][i+1] = x_cv+(x_cv-Bins['x'][i])

    Bins['Q2'] = np.zeros(Bins['Q2_cv'].shape[0]+1)
    Bins["Q2"][0] = Bins['Q2_cv'][0]-(Bins['Q2_cv'][1]-Bins['Q2_cv'][0])/2.
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

def GetCrossSections(Bins,sys_path="../expdata/data/xQ2binTable-xiaoxuan-060220+syst.npy"):

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
                else:
                    hist_CrossSection['non_empty'][key][iQ2*Bins['Nx']+iX] = True
                    hist_CrossSection['statunc'][key][iQ2*Bins['Nx']+iX] = np.sqrt(N)/lum_arg
                    hist_CrossSection['MCunc'][key][iQ2*Bins['Nx']+iX] = weight*Sample[key]['var_xsec']

    print("--Filling systematics in hist_CrossSection--")
    #! Barak's
    sys_bins = GetSysUnc(sys_path)
    sys_bins = sys_bins.loc[sys_bins['rs'] == rs]
    for key in hist_CrossSection['weights'].keys():
        for iQ2 in range(0, Bins['NQ2']):
            for iX in range(0, Bins['Nx']):
                if hist_CrossSection['Nevents'][key][iQ2*Bins['Nx']+iX] == 0:
                    hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = 0.
                else:
                    relative_sys = sys_bins.loc[sys_bins['x'] == Bins['x_cv'][iX]]
                    relative_sys = relative_sys.loc[sys_bins['Q2'] == Bins['Q2_cv'][iQ2]]
                    if relative_sys.iloc[0]['norm_c(%)']:
                        hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = np.sqrt((hist_CrossSection['xsecs'][key][iQ2*Bins['Nx']+iX]*relative_sys.iloc[0]['syst_u(%)']/100.)**2+(hist_CrossSection['xsecs'][key][iQ2*Bins['Nx']+iX]*relative_sys.iloc[0]['norm_c(%)']/100.)**2)
                    else:
                        hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = 0.

    """#! xiaoxuan
    #--- systematics uncertainties
    sys_bins = np.load(sys_path)
    for key in hist_CrossSection['weights'].keys():
        for iQ2 in range(0, Bins['NQ2']):
            for iX in range(0, Bins['Nx']):
                if hist_CrossSection['Nevents'][key][iQ2*Bins['Nx']+iX] == 0:
                    hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = 0.
                else:
                    relative_sys = np.where(
                        (sys_bins[:, 0] == Bins['x_cv'][iX]) & (sys_bins[:, 1] == Bins['Q2_cv'][iQ2]))
                    if np.any(relative_sys): #if bin has sys uncertainty
                        relative_sys = sys_bins[relative_sys][0, 2]
                    else: #if bin doesn't have sys uncertainty
                        relative_sys = 0

                    hist_CrossSection['sysunc'][key][iQ2*Bins['Nx']+iX] = hist_CrossSection['xsecs'][key][iQ2*Bins['Nx']+iX]*relative_sys
    """

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

def GetAnalysis(Bins, hist_CrossSection):
    hist_Analysis={}
    hist_Analysis['chi2s'] = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_Analysis['zscores'] = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_Analysis['x'] = np.zeros(Bins['Nx']*Bins['NQ2'])
    hist_Analysis['Q2'] = np.zeros(Bins['Nx']*Bins['NQ2'])

    hist_Analysis['xsecs'] = {}
    for key in hist_CrossSection['xsecs'].keys():
        hist_Analysis['xsecs'][key] = np.matrix(hist_CrossSection['xsecs'][key][hist_CrossSection['non_empty']['total']])

    print("--Filling hist_Analysis--")
    for iQ2 in range(0, Bins['NQ2']):
        for iX in range(0, Bins['Nx']):
            if hist_CrossSection['non_empty']['total'][iQ2*Bins['Nx']+iX]:
                ind = hist_CrossSection['non_empty']['map'][iQ2*Bins['Nx']+iX]

                max_d = hist_Analysis['xsecs']['max'][0,ind]
                min_d = hist_Analysis['xsecs']['min'][0,ind]

                hist_Analysis['chi2s'][iQ2*Bins['Nx']+iX] = (min_d-max_d)*inv_covmat_xsecs[ind,ind]*(min_d-max_d)
                hist_Analysis['zscores'][iQ2*Bins['Nx']+iX] = np.sqrt(hist_Analysis['chi2s'][iQ2*Bins['Nx']+iX])
                if hist_Analysis['zscores'][iQ2*Bins['Nx']+iX]<1:
                    hist_Analysis['zscores'][iQ2*Bins['Nx']+iX]=-1
                else:
                    hist_Analysis['zscores'][iQ2*Bins['Nx']+iX] = int(hist_Analysis['zscores'][iQ2*Bins['Nx']+iX])

            hist_Analysis['x'][iQ2*Bins['Nx']+iX] = Bins['x_cv'][iX]
            hist_Analysis['Q2'][iQ2*Bins['Nx']+iX] = Bins['Q2_cv'][iQ2]

            ##hist_Analysis['x'][iQ2*Bins['Nx']+iX] = Bins['x'][iX]+(Bins['x'][iX+1]-Bins['x'][iX])/2.
            ##hist_Analysis['Q2'][iQ2*Bins['Nx']+iX] = Bins['Q2'][iQ2]+(Bins['Q2'][iQ2+1]-Bins['Q2'][iQ2])/2.

    hist_Analysis['tot_chi2'] = (hist_Analysis['xsecs']['min']-hist_Analysis['xsecs']['max'])*inv_covmat_xsecs*(hist_Analysis['xsecs']['min']-hist_Analysis['xsecs']['max']).T
    hist_Analysis['tot_chi2'] /= Bins['total']

    hist_Analysis['zscores'][hist_Analysis['zscores'] >= 5] = 5

    #replace all 0 bins with nan for the hist to be properly colored:
    for key in hist_CrossSection['xsecs'].keys():
        hist_CrossSection['weights'][key][hist_CrossSection['weights'][key]==0] = np.nan
        hist_CrossSection['xsecs'][key][hist_CrossSection['xsecs'][key]==0] = np.nan
        hist_Analysis['chi2s'][hist_Analysis['chi2s']==0] = np.nan
        hist_Analysis['zscores'][hist_Analysis['zscores'] == 0] = np.nan

    hist_Analysis['zscores'][hist_Analysis['zscores'] == -1] = 0

    return hist_Analysis

def units(units):
    one=0.3893793721 #--GeV2 mbarn
    if   units=='GeV^-2':return 1
    elif units=='fb'    :return one*1e12 
    else: sys.exit('%s conversion not available')

if __name__ == "__main__":

    #--- Getting arguments from user
    if len(sys.argv) < 2:
        print "usage: ./chi2-test_generate.py <Luminosity (fb^-1)> <pdf_replica>"
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

    sub_sub_wdir = sub_wdir+"pdfrep" + str(rep)+"/"
    if not os.path.isdir(sub_sub_wdir): os.system("mkdir "+sub_sub_wdir)

    Sample_path = sub_sub_wdir+"Sample.p"
    hist_CrossSection_path = sub_sub_wdir+"hist_CrossSection.p"
    CovMat_path = sub_sub_wdir+"CovMat.p"
    hist_Analysis_path = sub_sub_wdir+"hist_Analysis.p"
    sys_path = pwd+"../expdata/src/ep_NC_optimistic-barak-100920.dat" #pwd+"../expdata/data/xQ2binTable-xiaoxuan-060220+syst.npy"

    Sample = {}
    hist_CrossSection = {}
    CovMat = {}
    hist_Analysis = {}
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
    rs = 140.71247279470288  # 140.7
    #-----------------------------------------------

    #--- acceptance/purity based binning
    Bins = GetBins(sys_path)

    #--- kinematic cuts
    def veto00(x,y,Q2,W2):
        if   W2 < 10       : return 0
        elif Q2 < 1        : return 0
        else               : return 1
        
    #--Getting events
    if not os.path.isfile(Sample_path):
        Sample = GetEvents(sub_sub_wdir, lum_arg, veto00, sign=1, rep=rep)
        save_obj(Sample,Sample_path)
    else:
        Sample = load_obj(Sample_path)

    #--Getting CrossSections
    if not os.path.isfile(hist_CrossSection_path):
        hist_CrossSection = GetCrossSections(Bins, sys_path=sys_path)
        save_obj(hist_CrossSection,hist_CrossSection_path)
    else:
        hist_CrossSection = load_obj(hist_CrossSection_path)

    #--Getting CovMat
    if not os.path.isfile(CovMat_path):
        CovMat = GetCovMat(Bins, hist_CrossSection)
        save_obj(CovMat,CovMat_path)
    else:
        CovMat = load_obj(CovMat_path)

    inv_covmat_xsecs = np.linalg.inv(CovMat['xsecs']['total'])

    #--Getting Analysis: Chi2s and Z-scores
    if not os.path.isfile(hist_Analysis_path):
        hist_Analysis = GetAnalysis(Bins, hist_CrossSection)
        save_obj(hist_Analysis,hist_Analysis_path)
    else:
        hist_Analysis = load_obj(hist_Analysis_path)

    print("All data has been generated here: "+sub_sub_wdir)
