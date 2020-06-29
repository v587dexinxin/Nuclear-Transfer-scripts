# -*- coding: utf-8 -*-
"""
Created on Fri Mar 01 20:53:15 2019

@author: xxli
"""

from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
cells = ['CCS','NT5','NT6','fESC']

tad_ccs = np.loadtxt(os.path.join(TadFolder , 'CCS_40K_allreps_union_small_domain.txt') , dtype = tad_type)



def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()
    

##-----------------------------------------CCS_standard--------------------------------------------------------------------------------
tad_new = []
for i in tad_ccs:
    chro = i['chr']
    if (i['end'] - i['start']) / 40000 >= 25:
        center = (i['start'] + i['end']) / 2
        s = center - 500000
        e = center + 500000
        tad_new.append((chro , int(s) , int(e)))        
Write2fils_nochr('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\CCS_standard_domain.txt' , tad_new)

##----------------------------------------Selected_Distance>=200K_domian-------------------------------------------------------------------

for c in cells:
    tads = np.loadtxt(os.path.join(TadFolder , c + '_40K_allreps.txt') , dtype = tad_type)
    tad_new = []
    for i in tads:
        chro = i['chr']
        if (i['end'] - i['start']) / 40000 >= 5:
            tad_new.append(i)
    Write2fils_nochr('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\' + c + '_domain_200K+.txt' , tad_new)
    
        
        
        
        
        

        