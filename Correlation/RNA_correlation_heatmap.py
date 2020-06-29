# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 15:18:21 2018

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Our Own Color Map
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')
                
                
RNAFolder = 'H:\\Workspace_New\\data\\RNA\\gene_expression'

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
    gene_all[i] = [0 , 0 , 0 , 0]
    for c in cell:
        ce = c.split("_")[0]
        mask = ([gene[c]['gene_name'] == i])
        overlap = gene[c][mask]
        if overlap.size == 0:
            fpkm = 0
        else:
            fpkm = overlap[0]['FPKM']
        gene_all[i][cells[ce]] += fpkm
    gene_all[i][0] = gene_all[i][0] / 3
    gene_all[i][1] = gene_all[i][1] / 4
    gene_all[i][2] = gene_all[i][2] / 3
    gene_all[i][3] = gene_all[i][3] / 3
    
    
    
RNA_all = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}            
for i in gene_name:
    if (gene_all[i][0] != 0) or (gene_all[i][1] != 0) or (gene_all[i][2] != 0) or (gene_all[i][3] != 0):
        if (gene_all[i][0] <= 3000) and (gene_all[i][1] <= 3000) and (gene_all[i][2] <= 3000) and (gene_all[i][3] <= 3000):
            RNA_all['CCS'].append(gene_all[i][0])
            RNA_all['NT5'].append(gene_all[i][1])
            RNA_all['NT6'].append(gene_all[i][2])
            RNA_all['fESC'].append(gene_all[i][3])
    
for k , v in RNA_all.items():
    v = np.array(v)
    RNA_all[k] = v
    
  
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
cor_matrix = np.zeros((4,4))
for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(RNA_all[i] , RNA_all[j])[0][1] , 3)
        
o = open('H:\\Workspace_New\\data\\correlation\\RNA_correlation_matrix_FPKM.txt' , 'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()

pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S3_figs\\S3A_RNA_FPKM_cor_100K.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=1,vmin = 0.91,cmap=my_cmap , origin = 'lower')

x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)
for i in range(len(x)):
    ax.plot([-0.5 , 3.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 3.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 20)


ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()


##hierarchy_cluster
data = np.vstack((RNA_all['CCS'] , RNA_all['NT5'] , RNA_all['NT6'] , RNA_all['fESC']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)


