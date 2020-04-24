#!/usr/bin bash

lhapdf=<path-to-lhapdf-install>
export PATH=$lhapdf/bin:$PATH
export PYTHONPATH=$PYTHONPATH:$lhapdf/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=$lhapdf/lib

path=`pwd`
export PYTHONPATH=$path/theory:$PYTHONPATH
export LHAPDF_DATA_PATH=$path/stf-grids


