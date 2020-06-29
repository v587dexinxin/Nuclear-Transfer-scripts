# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 15:38:21 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

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
                                            ['lightcyan' , '#FFFFFF' , 'lavender'])
my_cmap.set_bad('#D3D3D3')
                
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18 , 'X':19}
PCFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new'
#OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})

PC_all = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}

for c in cell:
    PCFil = c + '_compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCData = np.loadtxt(PCSource , dtype = sig_type)
    PC_all[c] = PCData['score']
            
    

    
cor_matrix = np.zeros((4,4))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(PC_all[i] , PC_all[j])[0][1] , 3)
        
o = open('D:\\Workspace_New\\data\\correlation\\Insulation_score_correlation_matrix.txt' , 'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()

f = open('D:\\Workspace_New\\data\\correlation\\RNA_correlation_matrix_log2(FPKM).txt' , 'r')
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

pp = PdfPages('D:\\Workspace_New\\Plot\\Correlation-Heatmap\\RNA_correlation_gene_3_(log2FPKM).pdf')
left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=1,vmin = cor_matrix.min(),cmap = my_cmap , origin = 'lower')

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

#ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
#fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()