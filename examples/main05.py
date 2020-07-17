#!/usr/bin/env python
import sys,os
sys.path.append(os.path.dirname( os.path.dirname(os.path.abspath(__file__) ) ) )
import numpy as np

import lhapdf
import argparse
#--matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
import pylab  as py

#--index conventions:
#--900: F2g  (NC F2 gamma channel) 
#--901: FLg  (NC FL gamma channel)
#--902: F2gZ (NC F2 gamma/Z channel)
#--903: FLgZ (NC FL gamma/Z channel)
#--904: F3gZ (NC F3 gamma/Z channel)
#--905: F2Z  (NC F2 Z channel)
#--906: FLZ  (NC FL Z channel)
#--907: F3Z  (NC F3 Z channel)
#--908: F2   (NC F2) 
#--909: FL   (NC FL) 
#--910: F3   (NC F3) 
#--930: W2m  (CC F2^{W-}) 
#--931: WLm  (CC FL^{W-})
#--932: W3m  (CC F3^{W-})
#--940: W2p  (CC F2^{W+})
#--941: WLp  (CC FL^{W+})
#--942: W3p  (CC F3^{W+})


Q2=[10,1000,10000]
X1=10**np.linspace(-4,-1)
X2=np.linspace(0.101,0.99)
X=np.append(X1,X2)

def get_data(target):
    #--set up dictionary based on target
    data = {_:{} for _ in ['idx','color','label','nrep']}

    if target=='proton':

        tag = 'p'
        tabnames=[        
                   'CT18ptxg'
                  ,'NNPDF31_lo_as_0118_SF'
                  ,'NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'
                  ,'JAM4EIC_p'
                 ]
        
        data['idx']['CT18ptxg'] = [900,901,902,903,904,905,906,907,908,909,910,930,931,932]
        data['idx']['NNPDF31_lo_as_0118_SF'] = [908,909,910,940,941,942,930,931,932]
        data['idx']['NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'] = [900,901,902,903,904,905,906,907,908,909,910,940,941,942,930,931,932]
        data['idx']['JAM4EIC_p'] = [900,901,902,903,904,905,906,907,908,909,910,940,941,942,930,931,932]
       
        data['nrep']['CT18ptxg'] =  0
        data['nrep']['NNPDF31_lo_as_0118_SF'] = 100
        data['nrep']['NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'] = 0
        data['nrep']['JAM4EIC_p'] = 50
 
        data['color']['CT18ptxg'] =  'red'
        data['color']['NNPDF31_lo_as_0118_SF'] = 'blue'
        data['color']['NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'] = 'orange'
        data['color']['JAM4EIC_p'] = 'green'

        data['label']['CT18ptxg'] =  'CT18'
        data['label']['NNPDF31_lo_as_0118_SF'] = 'NNPDF31 LO'
        data['label']['NNPDF31_nnlo_pch_as_0118_NC_Wm_Wp_SF'] = 'NNPDF31 NNLO'
        data['label']['JAM4EIC_p'] = 'JAM'

    elif target=='deuteron':

        tag = 'd'
        tabnames=[        
                  'JAM4EIC_d'
                 ]
        
        data['idx']['JAM4EIC_d'] = [900,901,902,903,904,905,906,907,908,909,910]
        
        data['nrep']['JAM4EIC_d'] = 50

        data['color']['JAM4EIC_d'] = 'green'

        data['label']['JAM4EIC_d'] = 'JAM'

    elif target=='helium':

        tag = 'h'
        tabnames=[        
                  'JAM4EIC_h'
                 ]
        
        data['idx']['JAM4EIC_h'] = [900,901,902,903,904,905,906,907,908,909,910]
        
        data['nrep']['JAM4EIC_h'] = 50

        data['color']['JAM4EIC_h'] = 'green'
        data['label']['JAM4EIC_h'] = 'JAM'

    else:
        print('Target: %s is not available!'%target)

    return data,tabnames,tag

def plot_stf(target):
   
    data,tabnames,tag = get_data(target) 
  
    #--collect data from different groups 
    for tabname in tabnames:
        data[tabname] = {}
        nrep = data['nrep'][tabname]
        idxs = data['idx'][tabname]
        for q2 in Q2: data[tabname][q2] = {_:[] for _ in idxs}
        for i in range(nrep+1): 
            stf=lhapdf.mkPDF(tabname,i)
            for q2 in Q2:
                for idx in idxs:
                    F=np.array([x*stf.xfxQ2(idx,x,q2) for x in X])
                    data[tabname][q2][idx].append(F)


    #--collect data from different groups 
    #for tabname in tabnames:
    #    data[tabname] = {}
    #    nrep = data['nrep'][tabname]
    #    idxs = data['idx'][tabname]
    #    for q2 in Q2:
    #        data[tabname][q2] = {_:[] for _ in idxs}
    #        for i in range(nrep+1): 
    #            stf=lhapdf.mkPDF(tabname,i)
    #            for idx in idxs:
    #                F=np.array([x*stf.xfxQ2(idx,x,q2) for x in X])
    #                data[tabname][q2][idx].append(F)
  
    #--plot data
 
    #--all channels 
    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 908 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][908][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 909 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][909][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,0.3)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 910 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F3=data[tabname][q2][910][i]
                ax.plot(X,F3,color=data['color'][tabname],alpha=alpha)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        #ax.set_ylim(0,0.9)
        if cnt==3: ax.text(0.1,0.7,r'$xF_3^{(%s)}$'%tag,size=35,transform=ax.transAxes)
        if cnt==9: ax.set_xlabel(r'$x$',size=20)
    
    py.tight_layout()
    filename = 'gallery/%s/FX.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()
       
    
    
    #--gamma channel
    nrows,ncols=3,2
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 900 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][900][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)\gamma}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 901 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][901][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,0.3)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)\gamma}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
    
    py.tight_layout()
    filename = 'gallery/%s/FXg.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()
    
    #--Z channel
    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 905 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][905][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 906 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][906][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,0.3)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 907 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F3=data[tabname][q2][907][i]
                ax.plot(X,F3,color=data['color'][tabname],alpha=alpha)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        #ax.set_ylim(0,0.9)
        if cnt==3: ax.text(0.1,0.7,r'$xF_3^{(%s)Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==9: ax.set_xlabel(r'$x$',size=20)
    
    
    py.tight_layout()
    filename = 'gallery/%s/FXZ.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()
    
    
    
    
    #--gamma/Z channel
    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 902 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][902][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)\gamma Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 903 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][903][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,0.3)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)\gamma Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 904 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F3=data[tabname][q2][904][i]
                ax.plot(X,F3,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,0.25)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==3: ax.text(0.1,0.7,r'$xF_3^{(%s)\gamma Z}$'%tag,size=35,transform=ax.transAxes)
        if cnt==9: ax.set_xlabel(r'$x$',size=20)
    
    
    py.tight_layout()
    filename = 'gallery/%s/FXgZ.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()
   
    #--no charged current for deuteron/helium
    if target=='deuteron': return
    if target=='helium':   return

 
    #--W plus
    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 940 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][940][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,12)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)W^+}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 941 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][941][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2.5)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)W^+}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 942 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F3=data[tabname][q2][942][i]
                ax.plot(X,F3,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,20)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==3: ax.text(0.1,0.7,r'$xF_3^{(%s)W^+}$'%tag,size=35,transform=ax.transAxes)
        if cnt==9: ax.set_xlabel(r'$x$',size=20)
    
    
    py.tight_layout()
    filename = 'gallery/%s/WXp.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()
    
    
    #--W minus
    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*3))
    
    cnt=0
    for q2 in Q2:
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 930 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F2=data[tabname][q2][930][i]
                ax.plot(X,F2,label=label,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,12)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==1: ax.legend(loc='lower left',bbox_to_anchor=(0,0.1),frameon=False,fontsize=15)
        if cnt==1: ax.text(0.1,0.7,r'$xF_2^{(%s)W^-}$'%tag,size=35,transform=ax.transAxes)
        if cnt==7: ax.set_xlabel(r'$x$',size=20)
        ax.set_ylabel(r'$Q^2=%0.1f~{\rm GeV^2}$'%q2,size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 931 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                FL=data[tabname][q2][931][i]
                ax.plot(X,FL,color=data['color'][tabname],alpha=alpha)
        #ax.set_ylim(0,2.5)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        if cnt==2: ax.text(0.1,0.7,r'$xF_L^{(%s)W^-}$'%tag,size=35,transform=ax.transAxes)
        if cnt==8: ax.set_xlabel(r'$x$',size=20)
    
        cnt+=1
        ax=py.subplot(nrows,ncols,cnt)
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        for tabname in tabnames:
            if 932 not in data[tabname][q2]: continue
            for i in range(data['nrep'][tabname]+1):
                if i == 0: label, alpha = data['label'][tabname], 1.0
                else:      label, alpha = None, 0.1
                F3=data[tabname][q2][932][i]
                ax.plot(X,F3,color=data['color'][tabname],alpha=alpha)
        ax.semilogx()
        ax.set_xticks([10e-5,10e-4,10e-3,10e-2,10e-1,1])
        ax.axhline(0,0,1,ls='--',color='black',alpha=0.5)
        #ax.set_ylim(0,20)
        if cnt==3: ax.text(0.1,0.7,r'$xF_3^{(%s)W^-}$'%tag,size=35,transform=ax.transAxes)
        if cnt==9: ax.set_xlabel(r'$x$',size=20)
    
    
    py.tight_layout()
    filename = 'gallery/%s/WXm.png'%target
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.clf()

if __name__=="__main__":

    ap = argparse.ArgumentParser()
    #--available targets: proton, deuteron, helium
    ap.add_argument('-t','--target',type=str, default='proton')

    args = ap.parse_args()
    plot_stf(args.target) 









