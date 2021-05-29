# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 20:05:00 2021

@author: xxli
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

                

##-------------------------------------------heatmap_new---------------------------------------------

cells = {'CCS':0 , 'MEF':1 , 'NTs':2 , 'fESC':3 , 'IPSC_P3':4 , 'IPSC_P20':5 , 'E14':6}
cell = ['CCS' , 'MEF' , 'NTs' , 'fESC' , 'IPSC_P3' , 'IPSC_P20' , 'E14']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18}
PCFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new'
#OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U4',np.float]})

PC_all = {}

for c in cell:
    PCFil = c + '_Traditonal_PC_200K_Compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCData = np.loadtxt(PCSource , dtype = sig_type)
    PC_all[c] = []
    for g in chrom:
        pc_data = PCData[PCData['chr'] == g]
        
        PC_all[c].extend(list(pc_data['score']))
            
    

    
cor_matrix = np.zeros((7,7))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(PC_all[i] , PC_all[j])[0][1] , 3)
        
        
        
pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\correlation\\Compartment_Correlation_cellline_heatmap_200K_7cellines.pdf')
left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
vmax = 0.99
vmin = 0.7
im = ax.imshow(cor_matrix , vmax=vmax , vmin = vmin , cmap = my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of Compartment between different cells',fontsize=25)

for i in range(len(x)):
    ax.plot([-0.5 , 6.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 6.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 20)

##Colorbar
ax = fig.add_axes([left + width + 0.035 , HB , 0.035 , 0.1])
cbar = fig.colorbar(im,cax = ax, orientation='vertical')
cbar.set_ticks([vmin , vmax])

pp.savefig(fig)
pp.close()

##hierarchy_cluster
data = np.vstack((PC_all['CCS'] , PC_all['NT5'] , PC_all['NT6'] , PC_all['F40'] , PC_all['F35']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)

