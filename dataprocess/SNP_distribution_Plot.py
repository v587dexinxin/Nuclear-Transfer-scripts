# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 22:01:17 2019

@author: han-luo
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


def Sig_To_200K(fil):
    """
    """
    data_type = np.dtype({'names' : ['chr' , 'pos'] , 
                      'formats' : ['S64' , np.int]})
    Data = np.loadtxt(fil , dtype = data_type , usecols = (1 , 2) , skiprows = 1)
    
    New_Data = {}
    for g in set(Data['chr']):
        New_Data[g] = {}
        tmp_data = Data[Data['chr'] == g]
        max_ = tmp_data['pos'].max()
        bin_size = max_ // 200000 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            data = line['pos'] // 200000
            if data >= bin_size:
                continue
            New_Data[g][data] += 1
        
    return New_Data
    
    
dataFil = 'C:\\Users\\lenovo\\Desktop\\PICTURE_New\\SNP\\GM_Snps.bed'

Sigdata = Sig_To_200K(dataFil)


sig = []
tricks = [0]
n = 0
for g in [str(x) for x in range(1,23)]:
    for i in Sigdata[g]:
        n += 1
        sig.append(i)
    tricks.append(n)
sig_index = np.arange(len(sig))
sig = np.array(sig)     


fig = plt.figure(figsize = (10 , 4))
ax = fig.add_axes([0.1  , 0.3 , 0.8 , 0.3])
ax.fill_between(sig_index , sig , facecolor = 'blue' , edgecolor = 'none' )
ax.set_ylim((-10 , 1000))
ax.set_xlim((0 , len(sig) + 200))
ax.set_xticks(tricks)
ax.set_xticklabels([str(x) for x in range(1,23)])




