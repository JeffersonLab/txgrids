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
from theory.mceg  import MCEG

#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

from scipy import stats
import pandas as pd
from collections import OrderedDict, Counter

wdir='.stats-tests'

#--physical params and cuts
rs= 140.7
lum='10:fb-1'
sign=1 #--electron=1 positron=-1

def veto00(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1


#--lhapdf set and stf idx
s = {'min':{}, 'max':{}}

s['min']['tabname'] = 'JAM4EIC'
s['min']['iset'] = 0
s['min']['iF2'] = 90001
s['min']['iFL'] = 90002
s['min']['iF3'] = 90003

s['max']['tabname'] = 'NNPDF31sx_nnlonllx_as_0118_LHCb_nf_6_SF'
s['max']['iset'] = 0
s['max']['iF2'] = 1001
s['max']['iFL'] = 1002
s['max']['iF3'] = 1003

#common keys
for key in s.keys():
    s[key]['wdir']      = wdir
    s[key]['sign']      = sign
    s[key]['rs'] = rs
    s[key]['fname'] = 'mceg00'
    s[key]['veto'] = veto00
    s[key]['fdata'] = wdir+"/"+s[key]['tabname']+"_data.po"


#--generate events
#output keys: ['Y', 'X', 'Q2', 'W', 'rs']
for key in s.keys():
    if not os.path.isfile(s[key]['fdata']):
        mceg=MCEG(**(s[key]))
        mceg.buil_mceg()
        ntot=10000
        s[key].update(mceg.gen_events(ntot))
        save(s[key], s[key]['fdata'])
    else:
        s[key] = load(s[key]['fdata'])

np.set_printoptions(threshold=sys.maxsize)

smin=np.array([s['min']['W'],s['min']['X'],s['min']['Q2']])
smax=np.array([s['max']['W'],s['max']['X'],s['max']['Q2']])

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
#--------------------------------------------------------------------------------------------------------------------






