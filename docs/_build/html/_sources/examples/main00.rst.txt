main00.py
=========

Gloals
------

- compute total cross section in :math:`ep` reaction

Code
----

.. code-block:: python

   #!/usr/bin/env python

   import sys,os
   sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
   import numpy as np
   from theory.tools import save, load
   from theory.idis  import IDIS
   
   def get_tot_xsec(tabname,rs=140.7,Q2min=1.0,W2min=10.0,neval=100000):
   
       def veto(x,y,Q2,W2):
           if   W2 < W2min  : return 0
           elif Q2 < Q2min  : return 0
           else             : return 1
   
       data={}
       data['tabname'] = tabname
       data['iset']    = 0
       data['iF2']     = 908
       data['iFL']     = 909    
       data['iF3']     = 910    
       data['sign']    = 1 #--electron=1 positron=-1
       data['veto']    = veto
       
       idis=IDIS(**data)
       
       data['neval'] = neval
       data['rs']    = rs
       data['iw']    = 0
       data['units'] = 'fb'
       data['mode']  = 'tot'
       
       val,err,Q = idis.get_cross_section(**data)
   
       print(' ')
       print('===============')
       print(tabname)
       print('%0.4e +/- %0.4e (%s) (Q=%f)'%(val,err,data['units'],Q))
       print(' ')
   
   
   if __name__=="__main__":
   
       get_tot_xsec('JAM4EIC',rs=140.7,Q2min=1.0,W2min=10.0)    
       get_tot_xsec('NNPDF31_lo_as_0118_SF',rs=140.7,Q2min=1.0,W2min=10.0)    
       get_tot_xsec('NNPDF31_nnlo_pch_as_0118_SF',rs=140.7,Q2min=1.0,W2min=10.0)    






