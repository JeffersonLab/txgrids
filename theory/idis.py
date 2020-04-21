#!/usr/bin/env python
import os,sys
import numpy as np
from   subprocess import Popen, PIPE
import pickle

import lhapdf
import vegas

import params as par
from   tools import checkdir,load,save,lprint

class IDIS:
    """
    Retrieves rows pertaining to the given keys from the Table instance
    represented by big_table.  Silly things may happen if
    other_silly_variable is not None.

    Args:
       big_table: An open Bigtable Table instance.
       keys: A sequence of strings representing the key of each table row
           to fetch.
       other_silly_variable: Another optional variable, that has a much
           longer name than the other args, and which does nothing.



    :param arg1: dictionary with parameters for initial setup
   
    .. note::

      An example of intersphinx is this: 
      you **cannot** use :mod:`pickle` on this class.
 
    
    Returns:
        A dict mapping keys to the corresponding table row data
        fetched. Each row is represented as a tuple of strings. For
        example:

        {'Serak': ('Rigel VII', 'Preparer'),
         'Zim': ('Irk', 'Invader'),
         'Lrrr': ('Omicron Persei 8', 'Emperor')}

        If a key from the keys argument is missing from the dictionary,
        then that row was not found in the table.

    :example:
    data={}
    data['tabname'] = tabname
    data['iset']    = iset   
    data['iF2']     = iF2    
    data['iFL']     = iFL    
    data['iF3']     = iF3    
    data['sign']    = 1 #--electron=1 positron=-1
    data['veto']    = veto

    
    Parameters
    ----------
    arg1 : int
        Description of arg1
    arg2 : str
        Description of arg2



    """

    def __init__(self,tabname,iset,iF2,iFL,iF3,sign,veto):

        self.tabname  =  tabname
        self.iset     =  iset   
        self.iF2      =  iF2    
        self.iFL      =  iFL    
        self.iF3      =  iF3    
        self.sign     =  sign   
        self.veto     =  veto   

        self.stf=lhapdf.mkPDF(self.tabname,self.iset)
        self.integ = vegas.Integrator([[0,1],[0,1]])

        if self.veto==None: self.veto=lambda x,y,Q2,W2:1

    def get_sigma_dxdy(self,x,y,iw):
        """
        positron sign =-1
        electron sign = 1
        """

        Q2=x*y*(self.rs**2-par.M2)
        F2=self.stf.xfxQ2(self.iF2,x,Q2) 
        FL=self.stf.xfxQ2(self.iFL,x,Q2) 
        F3=self.stf.xfxQ2(self.iF3,x,Q2) 

        if   iw==0: factor=1
        elif iw==1: factor=x
        elif iw==2: factor=y

        F1=(FL-(1+4*par.M2/Q2)*F2)/2/x
        factor*=4*np.pi*par.alfa**2/x/y/Q2
        return factor*((1-y-x**2*y**2*par.M2/Q2)*F2 + y**2*x*F1 +self.sign*(y-y**2/2)*x*F3)

    def get_sigma_dxdQ2(self,x,Q2,iw):
        """
        positron sign =-1
        electron sign = 1
        """
        y=Q2/x/(self.rs**2-par.M2)
        F2=self.stf.xfxQ2(self.iF2,x,Q2) 
        FL=self.stf.xfxQ2(self.iFL,x,Q2) 
        F3=self.stf.xfxQ2(self.iF3,x,Q2) 

        if   iw==0: factor=1
        elif iw==1: factor=x
        elif iw==2: factor=Q2

        F1=(FL-(1+4*par.M2/Q2)*F2)/2/x
        factor*=4*np.pi*par.alfa**2/x/y/Q2/(self.rs**2-par.M2)/x
        return factor*((1-y-x**2*y**2*par.M2/Q2)*F2 + y**2*x*F1 +self.sign*(y-y**2/2)*x*F3)

    def tgrand_dxdy(self,X):
        jx    = self.xmax-self.xmin
        x     = self.xmin+X[0]*jx

        if self.iymin==False: self.ymin=self.Q2min/(self.rs**2-par.M2)/x
        jy    = self.ymax-self.ymin
        y     = self.ymin+X[1]*jy

        Q2    = x*y*(self.rs**2-par.M2)
        W2    = par.M2+Q2/x*(1-x)
        xsec  = self.get_sigma_dxdy(x,y,self.iw)
        jac   = jx*jy
        wgt   = self.veto(x,y,Q2,W2)
        return xsec*wgt*jac

    def tgrand_dxdQ2(self,X):
        jx    = self.xmax-self.xmin
        x     = self.xmin+X[0]*jx

        if self.iQ2max==False: self.Q2max=(self.rs**2-par.M2)*x
        jQ2   = self.Q2max-self.Q2min
        Q2    = self.Q2min+X[1]*jQ2

        y     = Q2/(self.rs**2-par.M2)/x
        W2    = par.M2+Q2/x*(1-x)
        xsec  = self.get_sigma_dxdQ2(x,Q2,self.iw)
        jac   = jx*jQ2
        wgt   = self.veto(x,y,Q2,W2)
        return xsec*wgt*jac

    def units(self,units):
        one=0.3893793721 #--GeV2 mbarn
        if   units=='GeV^-2':return 1
        elif units=='fb'    :return one*1e12 
        else: sys.exit('%s conversion not available')

    def get_cross_section(self,neval,rs,mode,iw,
        xmin=None,xmax=None,ymin=None,ymax=None,iymin=None,
        Q2min=None,Q2max=None,iQ2min=None,units=None,**kargs):

        self.rs = rs 
        self.iw = iw 

        if mode=='xy':
      
            self.xmin  = xmin
            self.xmax  = xmax
            self.ymin  = ymin
            self.ymax  = ymax
            self.iymin = True
            result = self.integ(self.tgrand_dxdy, nitn=10, neval=neval)

        elif mode=='xQ2':
      
            self.xmin  = xmin
            self.xmax  = xmax
            self.Q2min = Q2min
            self.Q2max = Q2max
            self.iQ2max= True
            result = self.integ(self.tgrand_dxdQ2, nitn=10, neval=neval)

        elif mode=='tot':

            self.Q2min = 1 
            self.xmax  = 0.9
            self.ymax  = 1
            self.xmin  = self.Q2min/(self.rs**2-par.M2)
            self.iQ2max = False
            self.iymin  = False
            result = self.integ(self.tgrand_dxdy, nitn=10, neval=neval)
            #result = self.integ(self.tgrand_dxdQ2, nitn=10, neval=neval)

        return result.val*self.units(units)

if __name__=="__main__":

    
    tabname='JAM4EIC'             
    iset,iF2,iFL,iF3=0,90001,90002,90003  
    
    def veto(x,y,Q2,W2):
        if   W2 < 10          : return 0
        elif Q2 < 1           : return 0
        else                  : return 1
    
    data={}
    data['tabname'] = tabname
    data['iset']    = iset   
    data['iF2']     = iF2    
    data['iFL']     = iFL    
    data['iF3']     = iF3    
    data['sign']    = 1 #--electron=1 positron=-1
    data['veto']    = veto
    
    idis=IDIS(**data)

    data['neval'] = 10000
    data['rs']    = 140.7
    data['iw']    = 0
    data['Q2min'] = 1
    data['units'] = 'fb'
    #data['units'] = 'GeV^-2'

    data['mode']  = 'xy'
    data['xmin']  = 0.01
    data['xmax']  = 0.02
    data['ymin']  = 0.7
    data['ymax']  = 0.8
    
    print idis.get_cross_section(**data)

    #elif mode=='xQ2':
    #
    #    self.xmin  = data['xmin']
    #    self.xmax  = data['xmax']
    #    self.Q2min = data['Q2min']
    #    self.Q2max = data['Q2max']
    #    result = self.integ(tgrand_dxdQ2, nitn=10, neval=neval)

    #elif mode=='tot':





