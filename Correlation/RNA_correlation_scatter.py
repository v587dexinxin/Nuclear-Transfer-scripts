# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 19:12:39 2018

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
                                            ['#FFFFFF' ,'#CD0000'])
my_cmap.set_bad('#D3D3D3')
                
                
RNAFolder = 'D:\\Workspace_New\\data\\RNA\\gene_expression'

cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS_R1' , 'CCS_R2' , 'CCS_R3' , 'NT5_R1' , 'NT5_R2' , 'NT5_R3' , 'NT5_R4' ,'NT6_R1' , 'NT6_R2' , 'NT6_R3' , 'fESC_R1' , 'fESC_R2' , 'fESC_R3']
gene_type = np.dtype({'names':['gene_name','FPKM'],
                      'formats':['S64',np.float]})


gene = {}
gene_name = []
for c in cell:
    RNAFil = c + '.txt'
    RNASource = os.path.join(RNAFolder , RNAFil)
    RNAData = np.loadtxt(RNASource , skiprows = 1 , usecols = (1 , 7) , dtype = gene_type)
    gene[c] = RNAData
    for i in RNAData:
        if i['gene_name'] != '-':
            gene_name.append(i['gene_name'])
        
gene_name = set(gene_name)

gene_all = {}
for i in gene_name:
    gene_all[i] = []
    for c in cell:
        mask = ([gene[c]['gene_name'] == i])
        overlap = gene[c][mask]
        if overlap.size == 0:
            fpkm = 0
        else:
            fpkm = overlap[0]['FPKM']
        gene_all[i].append(fpkm)

    
RNA_all = {'CCS_R1':[] , 'CCS_R2':[] , 'CCS_R3':[] , 'NT5_R1':[] , 'NT5_R2':[] , 'NT5_R3':[] , 'NT5_R4':[] , 
           'NT6_R1':[] , 'NT6_R2':[] , 'NT6_R3':[] , 'fESC_R1':[] , 'fESC_R2':[] , 'fESC_R3':[]}            
for i in gene_name:
    RNA_all['CCS_R1'].append(np.log2(gene_all[i][0] + 1))
    RNA_all['CCS_R2'].append(np.log2(gene_all[i][1] + 1))
    RNA_all['CCS_R3'].append(np.log2(gene_all[i][2] + 1))
    RNA_all['NT5_R1'].append(np.log2(gene_all[i][3] + 1))
    RNA_all['NT5_R2'].append(np.log2(gene_all[i][4] + 1))
    RNA_all['NT5_R3'].append(np.log2(gene_all[i][5] + 1))
    RNA_all['NT5_R4'].append(np.log2(gene_all[i][6] + 1))
    RNA_all['NT6_R1'].append(np.log2(gene_all[i][7] + 1))
    RNA_all['NT6_R2'].append(np.log2(gene_all[i][8] + 1))
    RNA_all['NT6_R3'].append(np.log2(gene_all[i][9] + 1))
    RNA_all['fESC_R1'].append(np.log2(gene_all[i][10] + 1))
    RNA_all['fESC_R2'].append(np.log2(gene_all[i][11] + 1))
    RNA_all['fESC_R3'].append(np.log2(gene_all[i][12] + 1))
    
for k , v in RNA_all.items():
    v = np.array(v)
    RNA_all[k] = v
    

all_cells = ['CCS_R1ttCCS_R2' , 'CCS_R1ttCCS_R3' , 'CCS_R2ttCCS_R3' , 'NT5_R1ttNT5_R2' , 'NT5_R1ttNT5_R3' , 'NT5_R1ttNT5_R4' , 
             'NT5_R2ttNT5_R3' , 'NT5_R2ttNT5_R4' , 'NT5_R3ttNT5_R4' , 'NT6_R1ttNT6_R2' , 'NT6_R1ttNT6_R3' , 'NT6_R2ttNT6_R3' , 
             'fESC_R1ttfESC_R2' , 'fESC_R1ttfESC_R3' , 'fESC_R2ttfESC_R3']

size = (12, 12)
Left = 0.1 ; HB = 0.1 ; width = 0.8 ; HH = 0.8

pp = PdfPages('D:\\Workspace_New\\Plot\\Correlation-Heatmap\\replicated_scatter\\RNA_log2(FPKM).pdf')
cor_data = {}
for c in all_cells:
    c1 = c.split("tt")[0]
    c2 = c.split("tt")[1]
    cor = round(np.corrcoef(RNA_all[c1] , RNA_all[c2])[0][1],5)
    cor_data[c] = cor
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(RNA_all[c1] , RNA_all[c2] , alpha = 0.6 , c = 'red')
    ax.set_xlabel(c1 , size = 23)
    ax.set_ylabel(c2 , size = 23)
    ax.set_xlim(0,14)
    ax.set_ylim(0,14)
    ax.text(1 , 13 , 'R = ' + str(cor) , size = 20 )
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
    pp.savefig(fig) 
pp.close()       
