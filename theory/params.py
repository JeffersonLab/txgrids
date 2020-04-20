from mpmath import fp
from scipy.special import gamma
import numpy as np


#--set_constants
CA       = 3.0
CF       = 4.0/3.0
TR       = 0.5
TF       = 0.5
euler    = fp.euler 

#--masses
me       = 0.000511
mmu      = 0.105658
mtau     = 1.77684
mu       = 0.055
md       = 0.055
ms       = 0.2
mc       = 1.28
mb       = 4.18
mZ       = 91.1876
mW       = 80.398
M        = 0.93891897
Mpi      = 0.13803
Mk       = 0.493677
Mdelta   = 1.232

me2      = me**2 
mmu2     = mmu**2 
mtau2    = mtau**2
mu2      = mu**2  
md2      = md**2  
ms2      = ms**2  
mc2      = mc**2  
mb2      = mb**2  
mZ2      = mZ**2  
mW2      = mW**2  
M2       = M**2  
Mpi2     = Mpi**2  
Mdelta2  = Mdelta**2

#--set electro-week couplings
c2w      = mW2/mZ2
s2w      = 1.0-c2w
s2wMZ    = 0.23116
alfa     = 1/137.036
GF       = 1.1663787e-5   # 1/GeV^2
  
#--set stron coupling
alphaSMZ = 0.118


