# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 19:18:05 2018

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


#---------------------------------------------------Functions-------------------------------------------------------------
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'y', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
def caxis_PCA(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0.9)
    ax.spines['left'].set_linestyle('dotted')

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5', '6' ,'7' ,'8' , '9' , '10' ,'11' , '12' ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19']
R = 100000

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})


ATACFolder = 'H:\\Workspace_New\\data\\ATAC\\signal\\normalization\\bedgraph_100bp'
Outfolder = 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S3_figs'
OutFil = 'S3A_ATAC_bedgraph_100bp_cor_100K.pdf'

f = open('G:\\data\\genome\\mm10.txt' , 'r')
mm = {}
for i in f:
    i = i.strip().split()
    mm[i[0]] = int(i[1])
f.close()
    
#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.05 ; HB = 0.25 ; width = 0.5 ; HH = 0.5

all_data = {'CCS' : [] , 'NT5' : [] , 'NT6' : [] ,'fESC' : []}

for c in cell:
    ATACFil = c + '_100bp.bedgraph'
    ATACSource = os.path.join(ATACFolder , ATACFil)
    ATACData = np.loadtxt(ATACSource , dtype=sig_type_1)
    
    for g in chrom:
        atac = ATACData[ATACData['chr'] == g]
        length = mm[g] // R + 1
        new = np.zeros(length)
        for i in atac:
            start = i['start'] // R
            end = i['end'] // R
            score = i['score']
            n = end - start
            new[start] += score
            for j in np.arange(n):
                new[start + j] += score
        for i in new:
            all_data[c].append(i)
            
            

cor_matrix = np.zeros((4,4))

for i in cell:
    for j in cell:
        cor_matrix[cells[i] , cells[j]] = round(np.corrcoef(all_data[i] , all_data[j])[0][1] , 3)
        
o = open('H:\\Workspace_New\\data\\correlation\\ATAC_correlation_matrix.txt' , 'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()

f = open('H:\\Workspace_New\\data\\correlation\\ATAC_correlation_matrix.txt' , 'r')
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

pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S3_figs\\S3A_ATAC_bg100bp_cor_100K.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=1,vmin = 0.94,cmap=my_cmap)

x = ['CCS','NT5','NT6','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of Loop Strength  between different cells',fontsize=25)
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
data = np.vstack((all_data['CCS'] , all_data['NT5'] , all_data['NT6'] , all_data['fESC']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)


o = open('H:\\Workspace_New\\data\\correlation\\ATAC_bed_100Kb.txt' , 'w')
o.writelines('\t'.join(['CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
for i in range(len(all_data['CCS'])):
    o.writelines('\t'.join([str(all_data['CCS'][i]) , str(all_data['NT5'][i]) , str(all_data['NT6'][i]) , str(all_data['fESC'][i])]) + '\n')
o.close()

