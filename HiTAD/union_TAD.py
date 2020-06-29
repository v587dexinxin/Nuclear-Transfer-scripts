# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 20:05:01 2019

@author: xxli
"""

from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()


TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']
cells = ['CCS','NT5','NT6','fESC']


Tads_all = []
for c in cells:
    TadFil = np.loadtxt(os.path.join(TadFolder , c + '_40K_allreps_union_small_domain.txt') , dtype = tad_type)
    for i in TadFil:
        Tads_all.append(i)


Tads_all.sort(key = lambda x:(x[0],x[1],x[2]))

Tads_all = np.array(Tads_all ,dtype = tad_type)

union_tads = []
for g in chrom:
    Tads = Tads_all[Tads_all['chr'] == g]
    Tads_number = len(Tads)
    print Tads_number
    for i in range(Tads_number-1 , -1 , -1):
        start = Tads[i]['start']
        end = Tads[i]['end']
        prev_ss = Tads[i-1]['start'] - 40000
        prev_se = Tads[i-1]['start'] + 40000
        if (start >= prev_ss) and (start <= prev_se):
            Tads = list(Tads)
            Tads.remove(Tads[i])
            Tads = np.array(Tads,dtype = tad_type)
        else:
            pass
    print len(Tads)
    union_tads.extend(list(Tads))
        
union_tads = np.array(union_tads , dtype = tad_type)

Write2fils_nochr('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\union_domain.txt' , union_tads)            
    