Gettin Started
==============

dependencies
------------

- Only Linux and OSX is supported

- We recommend to install anaconda (python2) which 
  comes will all the necessary libraries


clonning the fipack repo
------------------------

Go to https://github.com/JeffersonLab/fitpack
and fork the repository under your GitHub account 
and clone your forked fitpack repo 

.. code-block:: shell

  git clone  git@github.com:<your-username>/fitpack.git

To update your fork with recent changes at the fitpack upstream 
you need to do the following 

.. code-block:: shell

   git remote add upstream git@github.com:JeffersonLab/fitpack.git
   git fetch upstream

This is done only once. After that you can sync the fork using 

.. code-block:: shell

   git pull upstream master


The fipack repository contains a minimalize version 
of the full program. The program needs to be extended
by cloning additional repos. To do that open the file 
``sync`` and uncomment the repos you are interested 
in having for a particular application. 
We have separated themaccording to 

- **root**: repos that a located at the root of the fitpack 

- **obslib**: collection of obserables 

- **grids**: mellin table grids for the observables

After chosing which repos you want to be installed, 
run the script ``sync``

.. code-block:: shell

  ./sync

This will clone only the repos that you have specified 
in the ``sync`` file. If you want to pull all the sub
repos, just run the sync script again and all the repos
will be updated.  

In addition, we recomment to not track the subrepos as
``submodules``. To avoid seen in the ``status`` output 
the changes on the submodules, use the lines ``git restore``.


You will  still see the untracked files

.. code-block:: shell

   Untracked files:
     (use "git add <file>..." to include in what will be committed)
           .gitmodules
           database/
           grids/
           nuclib/
           obslib/

but that is the normal behviour. 

Typically scripts such as those under fitlib, qcdlib, tools 
are typically modified independently by the user for testing/developing purposes 
(like adding print statements, etc). In some situation, generic 
features are needed to be added in those scripts. When that happen 
you can request a `pull request` using the GitHub system. 
The  fitpack admin will review the changes and evaluate if such changes 
can be incorporated under the master fitpack repo.


setups and checks
-----------------

Some environmental variables  need to be set. For csh add the following
lines on your `$HOME/.cshrc` 

.. code-block:: shell

   setenv FITPACK <path2>/fitpack
   setenv PYTHONPATH ${FITPACK}:${PYTHONPATH}
   setenv PATH ${FITPACK}/bin:${PATH}

and something analogous for bash.
To test if the JAM ecosystem is working, go to  fitlib and run 

.. code-block:: shell

  ./resman.py

If it crashes most likely some missing packges to be install via pip are 
needed 

.. code-block:: shell

  pip install <missing-package>


Notes
-----

- Each subrepois an independent repo. Once you clone the entire ecosystem
  you can go to each submodule and do pull/push.

- Single fits or MC runs should be made under a directory like workspace

The workspace 
-------------

- Use the `fitpack/workspace` directory as a template to start a new analysis from 
  scratch.

- `fitpack/workspace` is just a template, hence do not pull request any 
   changes on it.

- Instead, create a copy of such folder elsewere  or ask the mananger to create 
  a repo under JeffersonLab/JAM team


Next steps
----------

Checkout the tutorials




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


















