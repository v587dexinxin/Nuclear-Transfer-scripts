# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:53:51 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import pandas as pd

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
                
data_type = np.dtype({'names':['chr','start','r1','r2'],
                       'formats':['S8' , np.float , np.float , np.float]})
gene_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                       'formats':['S32' , 'S8' , 'S4' , np.int , np.int , np.float , np.float , np.float , np.float]})
cell = ['CCS_NT5' , 'CCS_NT6' , 'CCS_fESC']
R = 40000
boundaryFolder = 'H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\20181203stable_600K_diperse_diff_same'
geneFolder = 'H:\\Workspace_New\\data\\RNA\\gene_expression'
gene_FPKM = {'CCS_NT5' : [] , 'CCS_NT6' : [] , 'CCS_fESC' : []}
for c in cell:
    c1 = c.split("_")[0] ; c2 = c.split("_")[1]            
    boundaryFil = c1 + '_' + c2 + '_linear0.5_same.txt'
    boundarySource = os.path.join(boundaryFolder , boundaryFil)
    boundaryData = np.loadtxt(boundarySource , dtype = data_type)
    geneFil = 'all_gene_expression.txt'
    geneSource = os.path.join(geneFolder , geneFil)
    geneData = np.loadtxt(geneSource , skiprows = 1 , dtype = gene_type)
    
    for i in boundaryData:
        chro = i['chr']
        site = int(i['start'] * 1000000)
        start = site - 2 * R
        end = site + 2 * R
        gene = geneData[geneData['chr'] == 'chr' + chro]
        mask = (start <= gene['end']) & (end >= gene['start'])
        overlap = gene[mask]
        if overlap.size != 0:
            for j in overlap:
                if (j[c1] > 1) or (j[c2] > 1):
                    gene_FPKM[c].append((site , j['gene_name'] , j[c1] , j[c2]))
        else:
            continue
        
    d_type = np.dtype({'names':['boundary' , 'gene_name' , c1 , c2],
                       'formats':[np.int , 'S64' , np.float , np.float]})
    
    gene_FPKM[c] = np.array(gene_FPKM[c] , dtype = d_type)
    

pp = PdfPages('D:\\Workspace_New\\Plot\\Figure2\\boundary_gene_expression_boxplot\\diff_boundary_gene_expression.pdf')
data = [gene_FPKM['CCS_NT5']['CCS'] , gene_FPKM['CCS_NT5']['NT5'] , 
        gene_FPKM['CCS_NT6']['CCS'] , gene_FPKM['CCS_NT6']['NT6'] , 
        gene_FPKM['CCS_fESC']['CCS'] , gene_FPKM['CCS_fESC']['fESC']]
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.boxplot(data , showmeans=True , showfliers=False , positions=[0.1,0.8,2.1,2.8,4.1,4.8])
ax.set_xticklabels(['CCS' , 'NT5' , 'CCS' , 'NT6' , 'CCS' , 'fESC'] , fontsize = 28)
ax.set_yticks([0 , 10 , 20 , 30])
ax.set_yticklabels(['0' , '10' , '20' , '30'] , fontsize = 15)
pp.savefig(fig)
pp.close()
  
