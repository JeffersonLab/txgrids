#!/usr/bin/env python
import os
import sys
import numpy as np
from subprocess import Popen, PIPE
import pickle

import lhapdf
import vegas

import params as par
from tools import checkdir, load, save, lprint
from tools import convert_lum


class MCEG:
    """
    :param wdir:    directory to store the results
    :param tabname: The lhapdf  table name see LHAPDF grids
    :param iset:    The table set  
    :param iF2:     flavor index for F2
    :param iFL:     flavor index for FL
    :param iF3:     flavor index for F3
    :param sign:    1 for electron 1 for positron
    :param rs:      center of mass energy
    :param fname:   file name for the generator
    :param veto:    user defined cuts
    """

    def __init__(self, wdir, tabname, iset, iF2, iFL, iF3, sign, rs, fname, veto, **kargs):

        self.wdir = wdir
        self.tabname = tabname
        self.iset = iset
        self.iF2 = iF2
        self.iFL = iFL
        self.iF3 = iF3
        self.sign = sign
        self.rs = rs
        self.fname = fname
        self.veto = veto

        self.stf = lhapdf.mkPDF(self.tabname, self.iset)

        self.Q2min = 1
        self.Q2max = self.rs**2-par.M2
        self.xmin = self.Q2min/(self.rs**2-par.M2)
        self.xmax = 1
        self.ymax = 1

        self.iw = 0  # --observable

        self.integ = vegas.Integrator([[0, 1], [0, 1]])

        if self.veto == None:
            self.veto = lambda x, y, Q2, W2: 1

    def _get_sigma_dxdy(self, x, y, iw):

        Q2 = x*y*(self.rs**2-par.M2)
        F2 = self.stf.xfxQ2(self.iF2, x, Q2)
        FL = self.stf.xfxQ2(self.iFL, x, Q2)
        F3 = self.stf.xfxQ2(self.iF3, x, Q2)

        if iw == 0:
            factor = 1
        elif iw == 1:
            factor = x
        elif iw == 2:
            factor = Q2

        F1 = ((1+4*par.M2/Q2)*F2-FL)/2/x
        factor *= 4*np.pi*par.alfa**2/x/y/Q2
        return factor*((1-y-x**2*y**2*par.M2/Q2)*F2 + y**2*x*F1 + self.sign*(y-y**2/2)*x*F3)

    def _f(self, X):
        jx = self.xmax-self.xmin
        x = self.xmin+X[0]*jx
        ymin = self.Q2min/(self.rs**2-par.M2)/x
        jy = self.ymax-ymin
        y = ymin+X[1]*jy
        Q2 = x*y*(self.rs**2-par.M2)
        W2 = par.M2+Q2/x*(1-x)

        Q2max = x*self.ymax*(self.rs**2-par.M2)
        if Q2 > Q2max:
            return 0

        xsec = self._get_sigma_dxdy(x, y, self.iw)
        jac = jx*jy
        wgt = self.veto(x, y, Q2, W2)
        return xsec*wgt*jac

    def _vegas_integrate(self, neval=1000, iw=0, nitn=10):
        self.iw = iw
        result = self.integ(self._f, nitn=nitn, neval=neval)
        return result

    def buil_mceg(self, neval=1000, nitn=10):
        """
        build the vegas integrator as MCEG
        """

        print('integrating cross section...')
        cnt = 0
        while 1:
            print('trial %d' % cnt)
            result = self._vegas_integrate(neval=neval, nitn=nitn)
            if result.Q < 0.2:
                neval *= 10
            else:
                break
            cnt += 1
            if cnt == 1:
                break
        print(result.summary())

        sigtot = result.val
        Q = result.Q

        if Q < 0.1:
            msg = '(bad) '
        else:
            msg = '(good)'

        print('\n')
        print('sig tot = %.2e GeV^-2' % sigtot)
        print('Q       = %.2f %s    ' % (Q, msg))

        mceg_name = '%s/%s.po' % (self.wdir, self.fname)
        print('saving mceg at %s' % mceg_name)
        checkdir(self.wdir)
        with open(mceg_name, 'wb') as f:
            pickle.dump(self.integ, f)

    def _load_mceg(self):
        file = open('%s/%s.po' % (self.wdir, self.fname), 'rb')
        self.integ = pickle.load(file)

    def get_xsectot(self):
        """
        :return: total number of events for a given 
                 of the current generator and input luminosity
        """
        self._load_mceg()
        tot = 0
        for X, wgt in self.integ.random():
            tot += wgt*self._f(X)
        return tot

    def get_ntot(self, lum):
        """
        :param lum:  str e.g.  '10:fb-1' 
        :return: total number of events for a given 
                 of the current generator and input luminosity
        """
        self._load_mceg()
        tot = 0
        for X, wgt in self.integ.random():
            tot += wgt*self._f(X)
        lum = convert_lum(lum)
        return int(lum*tot)

    def _get_batch_size(self):
        self._load_mceg()
        cnt = 0
        for X, wgt in self.integ.random():
            cnt += 1
        return cnt

    def gen_events(self, nevents):
        """
        Generate events

        :param nevents: number of desired events
        :return: dictionary with arrays for X,Y,Q2,W.
                 it also includes root s used to generate
                 the mceg in case is needed for cross checks
        """
        print('\nstart event generation')

        self._load_mceg()
        bsize = self._get_batch_size()
        nruns = nevents/bsize+1

        X, Y, W = [], [], []

        cnt = 0
        for i in range(nruns):
            _X, _Y, _W = [], [], []
            for V, wgt in self.integ.random():
                _w = wgt*self._f(V)
                if _w <= 0:
                    continue
                cnt += 1
                lprint('%d/%d' % (cnt, nevents))
                jx = self.xmax-self.xmin
                x = self.xmin+V[0]*jx
                ymin = self.Q2min/(self.rs**2-par.M2)/x
                jy = self.ymax-ymin
                y = ymin+V[1]*jy

                _X = np.append(_X, x)
                _Y = np.append(_Y, y)
                _W = np.append(_W, _w)
            _W /= np.sum(_W)
            X = np.concatenate((X, _X), axis=None)
            Y = np.concatenate((Y, _Y), axis=None)
            W = np.concatenate((W, _W), axis=None)
        W /= nruns
        X = X[:nevents]
        Y = Y[:nevents]
        W = W[:nevents]
        print('')

        Q2 = X*Y*(self.rs**2-par.M2)

        data = {}
        data['X'] = X
        data['Y'] = Y
        data['Q2'] = Q2
        data['W'] = W
        data['rs'] = self.rs

        return data


if __name__ == "__main__":

    wdir = 'tmp'

    #--physical params
    rs = 140.7
    lum = '10:fb-1'
    sign = 1  # --electron=1 positron=-1
    W2cut = 10

    #--lhapdf set and stf idx
    tabname = 'JAM4EIC'
    iset, iF2, iFL, iF3 = 0, 90001, 90002, 90003

    def veto00(x, y, Q2, W2):
        if W2 < 10:
            return 0
        elif Q2 < 1:
            return 0
        else:
            return 1

    fname, veto = 'mceg00', veto00

    data = {}
    data['wdir'] = wdir
    data['tabname'] = tabname
    data['iset'] = iset
    data['iF2'] = iF2
    data['iFL'] = iFL
    data['iF3'] = iF3
    data['sign'] = sign
    data['rs'] = rs
    data['fname'] = fname
    data['veto'] = veto

    mceg = MCEG(**data)
    mceg.buil_mceg()
    ntot = mceg.get_ntot(lum)
    ntot = 10000
    mceg.gen_events(ntot)
