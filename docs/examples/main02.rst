main02.py
=========

Gloals
------

- compute bin/total cross section in :math:`ep` reaction


Code
----

First we include some headers

.. code-block:: python

   import sys,os
   sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
   import numpy as np
   from theory.tools import save, load
   from theory.idis  import IDIS
   

``sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )`` 
will add paths to be able to load the theory module. 

Next we setup lhapdf structure functions

.. code-block:: python
   tabname='JAM4EIC'             
   iset,iF2,iFL,iF3=0,90001,90002,90003  


Next we set physical parameters of the reaction

.. code-block:: python
   
   rs   = 140.7
   lum  = '10:fb-1'
   sign = 1  #--electron=1 positron=-1

Select lhapdf tables and flavor flags for structure functions

.. code-block:: python
   
   tabname='JAM4EIC'             
   iset,iF2,iFL,iF3=0,90001,90002,90003  

Define a user defined global cuts

.. code-block:: python
   
   def veto(x,y,Q2,W2):
       if   W2 < 10          : return 0
       elif Q2 < 1           : return 0
       else                  : return 1
   
Create a dictionary with all the parameters
and initialize the ``IDIS`` class

.. code-block:: python
   
   data={}
   data['tabname'] = tabname
   data['iset']    = iset   
   data['iF2']     = iF2    
   data['iFL']     = iFL    
   data['iF3']     = iF3    
   data['sign']    = 1 #--electron=1 positron=-1
   data['veto']    = veto

   idis=IDIS(data)

Next we setup the kimeatics to compute the
cross section 
   
.. code-block:: python
   
   data['neval'] = 10000
   data['rs']    = 140.7
   data['iw']    = 0
   data['units'] = 'fb'  # or 'GeV^-2'
   
   data['mode']  = 'xy'
   data['xmin']  = 0.01
   data['xmax']  = 0.02
   data['ymin']  = 0.7
   data['ymax']  = 0.8
   
   print(idis.get_cross_section(data))

Here use ``mode=xy`` to integrate cross sections in :math:`x-y` space.  The
parameter ``iw=0,1,2`` is to integrate with  a factor of :math:`1,x,y`
respectively. ``neval`` controls the vegas samples. 

The ``mode`` parameter can be changed to compute 
the cross section in :math:`x-Q^2` space e.g

.. code-block:: python

   data['mode']  = 'xQ2'
   data['xmin']  = 0.01
   data['xmax']  = 0.02
   data['Q2min'] = 5.0
   data['Q2max'] = 10.0
   
   print('%0.3e'%idis.get_cross_section(data))

or to compute the toal cross section 
with cuts

.. code-block:: python
   
   data['mode']  = 'tot'
   print(idis.get_cross_section(data))




