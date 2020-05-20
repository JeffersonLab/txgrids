#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np
from theory import tools 
from tools  import checkdir

class CONVERTER:

    def __init__(self,name,iflavs):

        fnames=os.listdir(name)
        self.info   = [_ for _ in fnames if _.endswith('.info')][0]
        self.isets  = [_ for _ in fnames if _.endswith('.info')==False]
        self.name   = name
        self.iflavs = iflavs        
        self.mod_info()
        self.mod_isets()

    def savefile(self,name,L):
        L=['%s\n'%_ for _ in L]
        F=open(name,'w')
        F.writelines(L)
        F.close()

    def mod_info(self):
        L=open('%s/%s'%(self.name,self.info)).readlines()
        L=[_.strip() for _ in L]
        for i in range(len(L)):
            if L[i].startswith('Flavors'):
                L[i]='Flavors: ['
                for iflav in self.iflavs: L[i]+='%d,'%iflav
                L[i]=L[i].rstrip(',')+']'
        self.savefile('%s/%s'%(self.name,self.info),L)

    def RepresentsInt(self,s):
        try: 
            int(s)
            return True
        except ValueError:
            return False

    def mod_isets(self):
        l=''
        for iflav in self.iflavs: l+='%d '%iflav

        nflavs=len(self.iflavs)

        for iset in self.isets:
            L=open('%s/%s'%(self.name,iset)).readlines()
            L=[_.strip() for _ in L]
            for i in range(len(L)):
                row=L[i].split()
                nl=len(row)
                if nl==nflavs:
                    if self.RepresentsInt(row[0]):
                        L[i]=l
            self.savefile('%s/%s'%(self.name,iset),L)

if __name__=="__main__":

    iflavs=[908, 909, 910, 
            940, 941, 942,
            930, 931, 932,]

    name='NNPDF31_lo_as_0118_SF'
    CONVERTER(name,iflavs)

    name='NNPDF31_nnlo_pch_as_0118_SF'
    CONVERTER(name,iflavs)





