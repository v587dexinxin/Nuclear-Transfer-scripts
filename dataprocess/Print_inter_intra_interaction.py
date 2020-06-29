# -*- coding: utf-8 -*-
"""
Created on Mon May 04 16:59:30 2020

@author: han-luo
"""


from __future__ import division
import numpy as np



for c in ['CCS_R1' , 'CCS_R2' , 'NT5_R1' , 'NT5_R2' , 'NT6_R1' , 'NT6_R2' , 'fESC_R1' , 'fESC_R2']:
    inter = 0 ; intra = 0
    data = open('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/' + c.split('_')[0] + '_workspace/' + c + '/Filtered_Bed/' + c + '_Valid.bed')
    for i in data:
        i = i.split()
        if (i[1] == i[8]) & (abs(int(i[3]) - int(i[10])) >= 20000):
            intra += 1
        elif i[1] != i[8]:
            inter += 1
    print inter , intra
        
