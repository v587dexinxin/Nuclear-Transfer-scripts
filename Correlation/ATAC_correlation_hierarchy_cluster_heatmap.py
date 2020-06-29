# -*- coding: utf-8 -*-
"""
Created on Sat Nov 02 21:42:29 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

                
f = open('H:\\Workspace_New\\data\\correlation\\ATAC_correlation_matrix_log2(RPKM).txt' , 'r')
a = f.readlines()
cor_matrix = np.zeros((4,4))
cor_matrix[0 , 0] = a[1].split()[1]
cor_matrix[0 , 1] = a[1].split()[2]
cor_matrix[0 , 2] = a[1].split()[3]
cor_matrix[0 , 3] = a[1].split()[4]
cor_matrix[1 , 0] = a[2].split()[1]
cor_matrix[1 , 1] = a[2].split()[2]
cor_matrix[1 , 2] = a[2].split()[3]
cor_matrix[1 , 3] = a[2].split()[4]
cor_matrix[2 , 0] = a[3].split()[1]
cor_matrix[2 , 1] = a[3].split()[2]
cor_matrix[2 , 2] = a[3].split()[3]
cor_matrix[2 , 3] = a[3].split()[4]
cor_matrix[3 , 0] = a[4].split()[1]
cor_matrix[3 , 1] = a[4].split()[2]
cor_matrix[3 , 2] = a[4].split()[3]
cor_matrix[3 , 3] = a[4].split()[4]
       
pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\Compartment_Correlation_heatmap_200K.pdf')
left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
vmax = 0.99
vmin = 0.88
im = ax.imshow(cor_matrix , vmax=vmax , vmin = vmin , cmap = my_cmap , origin = 'lower')

x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 40)
ax.set_yticklabels(x,fontsize = 40)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)

for i in range(len(x)):
    ax.plot([-0.5 , 3.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 3.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 40)

##Colorbar
ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
cbar = fig.colorbar(im,cax = ax, orientation='vertical')
cbar.set_ticks([vmin , vmax])

pp.savefig(fig)
pp.close()

##hierarchy_cluster
data = np.vstack((PC_all['CCS'] , PC_all['NT5'] , PC_all['NT6'] , PC_all['fESC']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)

