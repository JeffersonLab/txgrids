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

fig = plt.figure()
ax = fig.add_subplot(111)
#plt.tight_layout(rect=[0, 0.08, 1., 0.95])
ax.title.set_text('')

n1 = 10000
mu1 = 0.
sig1 = 1.
rvs1 = stats.norm.rvs(size=n1, loc=mu1, scale=sig1)


(n, bins, patches) = ax.hist(rvs1, bins='fd',density=True, facecolor='blue', edgecolor='blue', histtype="step", alpha=1., label='min', ls='-', linewidth=1.5)


kde = stats.gaussian_kde(rvs1)

n = m = 100 

bins = np.linspace(min(rvs1), max(rvs1), n)

ax.plot(bins, kde.pdf(bins), color='blue', ls='--')

plt.savefig("test_kde.pdf")
plt.clf

