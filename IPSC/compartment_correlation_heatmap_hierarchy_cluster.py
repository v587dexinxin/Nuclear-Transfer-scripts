# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 20:06:08 2019

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

                
cells = {'MEF':0 , 'IPS_P3':1 , 'E14':2}
cell = ['MEF' , 'IPS_P3' , 'E14']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18}
PCFolder = 'H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new'
#OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})

PC_all = {}

for c in cell:
    PCFil = c + '_Compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCData = np.loadtxt(PCSource , dtype = sig_type)
    PC_all[c] = PCData[PCData['chr'] != 'X']['score']
            
    

    
cor_matrix = np.zeros((3,3))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(PC_all[i] , PC_all[j])[0][1] , 3)
        
pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPS_Compartment_Correlation_heatmap_200K.pdf')

size_axes = [Left, HB, width, HH]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
vmax = 0.99
vmin = 0.85
im = ax.imshow(cor_matrix , vmax=vmax , vmin = vmin , cmap = my_cmap , origin = 'lower')

x = ['MEF','IPS_P3','E14']
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 40)
ax.set_yticklabels(x,fontsize = 40)
ax.set_title('Correlation of Compartment between different cells',fontsize=25)

for i in range(len(x)):
    ax.plot([-0.5 , 2.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 2.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 20)

##Colorbar
ax = fig.add_axes([Left + width + 0.035 , HB , 0.035 , 0.1])
cbar = fig.colorbar(im,cax = ax, orientation='vertical')
cbar.set_ticks([vmin , vmax])

pp.savefig(fig)
pp.close()

##hierarchy_cluster
data = np.vstack((PC_all['MEF'] , PC_all['IPS_P3'] , PC_all['E14']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)

