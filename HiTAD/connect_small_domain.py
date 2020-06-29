# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 21:50:32 2019

@author: xxli
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os
import math

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

cells = ['CCS','NT5','NT6','fESC']

TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' ,'X']

for c in cells:
    TadFil = c + '_40K_allreps.txt'
    TadSource = os.path.join(TadFolder , TadFil)
    TadData = np.loadtxt(TadSource , usecols = (0 , 1 , 2) ,dtype = tad_type)
    out = open('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\' + c + '_40K_allreps_union_small_domain.txt' , 'w')
    for g in chrom:
        tmp = []
        tads = TadData[TadData['chr'] == g]
        for i in range(len(tads)-1):
            tmp.append(tads[i])
            if tads[i + 1][1] -  tads[i][2] >= 80000:
                print tads[i] ,tads[i + 1]
                tmp.append((g , tads[i][2] , tads[i + 1][1]))
            else:
                continue
        tmp = Sort(tmp , 0 , 1)
        for i in tmp:
            out.writelines(i[0] + '\t' +  str(i[1]) + '\t' + str(i[2]) + '\n')
    out.close()
        
            
    