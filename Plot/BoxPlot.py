# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 21:20:24 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from  scipy.stats import ttest_rel
import pandas as pd


f = np.loadtxt('C:\\Users\\Administrator\\Desktop\\snp_numbers.txt' , usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20) ,
               dtype = np.int)


data = {'CCS_based_NT5':f[1] , 'CCS_based_NT6':f[0] , 'CCS_based_fESC':f[2]}
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
x = ax.boxplot([data['CCS_based_NT5'] , data['CCS_based_NT6'] , data['CCS_based_fESC']], showmeans=True , showfliers=False,patch_artist=True)
colors = ['pink', 'lightblue', 'lightgreen']
for patch, color in zip(x['boxes'], colors):
    patch.set_facecolor(color)
ax.set_xticklabels(['CCS_based_NT5' , 'CCS_based_NT6' , 'CCS_based_fESC'] , fontsize = 20)

ax.set_yticks([0 , 10 , 20 , 30])
ax.set_yticklabels(['0' , '10' , '20' , '30'] , fontsize = 20)
ax.set_ylabel('Numbers of SNP' , fontsize = 30)
plt.savefig('D:\\Workspace_New\\Plot\\Snp_numbers.png')

NT5_NT6 = ttest_rel(data['CCS_based_NT5'], data['CCS_based_NT6'])
NT5_fESC = ttest_rel(data['CCS_based_NT5'], data['CCS_based_fESC'])
NT6_fESC = ttest_rel(data['CCS_based_NT6'], data['CCS_based_fESC'])

