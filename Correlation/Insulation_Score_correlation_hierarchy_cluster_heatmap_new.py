# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:15:29 2020

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

                
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'F40':3 , 'F35':4}
cell = ['CCS' , 'NT5' , 'NT6' , 'F40' , 'F35']

ISFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\Insulation_score'
#OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U4',np.float]})

IS_all = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[] , 'F35':[]}

for c in cell:
    ISFil = c + '_Insulation_score_40K.txt'
    ISSource = os.path.join(ISFolder , ISFil)
    ISData = np.loadtxt(ISSource , dtype = sig_type , skiprows = 1 , usecols = (0 , 2))
    IS_all[c] = ISData[ISData['chr'] != 'X']['score']
            
    

    
cor_matrix = np.zeros((5,5))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(IS_all[i] , IS_all[j])[0][1] , 3)
        
pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\correlation\\Insulation_score_Correlation_heatmap_40K.pdf')
left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
vmax = 1
vmin = 0.9
im = ax.imshow(cor_matrix , vmax=vmax , vmin = vmin , cmap = my_cmap , origin = 'lower')

x = ['CCS','NT5','NT6','fESC','F35']
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 40)
ax.set_yticklabels(x,fontsize = 40)
ax.set_title('Correlation of Insulation Score between different cells',fontsize=25)

for i in range(len(x)):
    ax.plot([-0.5 , 4.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 4.5] , color="black")
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
data = np.vstack((IS_all['CCS'] , IS_all['NT5'] , IS_all['NT6'] , IS_all['fESC'] , IS_all['F35']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)

