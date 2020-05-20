#!/usr/bin/env python
######################################
# authors: N.Sato (nsato@jlab.org) 
#
# last update: 04-20-20
#####################################
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory.tools import save, load
from theory.idis  import IDIS

#--matplotlib

tabname='JAM4EIC'             
iset,iF2,iFL,iF3=0,90001,90002,90003  

def veto(x,y,Q2,W2):
    if   W2 < 10          : return 0
    elif Q2 < 1           : return 0
    else                  : return 1

data={}
data['tabname'] = tabname
data['iset']    = iset   
data['iF2']     = 908
data['iFL']     = 909    
data['iF3']     = 910    
data['sign']    = 1 #--electron=1 positron=-1
data['veto']    = veto

idis=IDIS(**data)

data['neval'] = 10000
data['rs']    = 140.7
data['iw']    = 0
data['units'] = 'fb'
#data['units'] = 'GeV^-2'

data['mode']  = 'tot'
val,err,Q = idis.get_cross_section(**data)
print('%0.4e +/- %0.4e (%s) (Q=%f)'%(val,err,data['units'],Q))

data['mode']  = 'xy'
data['xmin']  = 0.01
data['xmax']  = 0.02
data['ymin']  = 0.7
data['ymax']  = 0.8
val,err,Q = idis.get_cross_section(**data)
print('%0.4e +/- %0.4e (%s) (Q=%f)'%(val,err,data['units'],Q))


data['mode']  = 'xQ2'
data['xmin']  = 0.01
data['xmax']  = 0.02
data['Q2min'] = 5.0
data['Q2max'] = 10.0
val,err,Q = idis.get_cross_section(**data)
print('%0.4e +/- %0.4e (%s) (Q=%f)'%(val,err,data['units'],Q))




