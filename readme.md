# txgrids

Theoretical cross section grids for EIC YR

## dependencies 

- python 2 or 3
- lhapdf

## setup

- use the script  ``install.py`` to create symlinks to your 
  lhapdf directory 

## usage

Users should create excecutable scripts  at the folder ``main``
e.g. ``main/mainXY.py``  for some available index XY

For instance to run ``main/main00/py`` do the following:

``cd main``

``./main00.py``

## What is available 

| observable     | reaction             | group | main      | comments | contact |
| :--:           | :--:                 | :--:  | :--:      | :--:     | :--:    |
| ``dsig/dxdQ2`` | ``nu + p -> nu' +X`` | NNPDF | main00.py | in-dev   | [rojo]  |


[rojo]:<mailto:j.rojo@vu.nl>
[nsato]:<mailto:nsato@jlab.org>

## TODO

### General 

- add scripts to generate lhapdf grids for cross sections
- add vanilla plotting scripts and examples in jupyter-notebook

### CC proton

- currently only  CC  for ``e+p -> v+X`` is available from NNPDF
  we need to get NC for ``e+p -> v+X`` and from other groups


### PVDIS proton 

- neet grids for ``F_2,3 (gamma/Z)``
- neet grids for ``F_2,3 (W-)``
- neet grids for ``g_1,5 (gamma/Z)``
- neet grids for ``g_1,5 (W-)``











