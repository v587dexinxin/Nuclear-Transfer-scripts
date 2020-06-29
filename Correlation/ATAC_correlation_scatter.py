# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 19:55:43 2018

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


# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

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
chrom=['1' , '2' , '3' , '4' , '5', '6' ,'7' ,'8' , '9' , '10' ,'11' , '12' ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X']
R = 5000

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})


ATACFolder = 'D:\\Workspace_New\\data\\ATAC\\signal\\normalization\\bedgraph_100bp'
Outfolder = 'D:\\Workspace_New\\Plot\\Correlation-Heatmap\\replicated_scatter'
OutFil = 'ATAC_5K.pdf'

f = open('G:\\data\\genome\\mm10.txt' , 'r')
mm = {}
for i in f:
    i = i.strip().split()
    mm[i[0]] = int(i[1])
f.close()
    
#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6


all_data = {'CCS' : [] , 'CCS1' : [] , 'CCS2':[] , 'NT5' : [] , 'NT51':[] , 'NT52':[] , 
            'NT6' : [] , 'NT61' : [] , 'NT62':[] , 'fESC' : [] , 'fESC1' : [] , 'fESC2' : []}


for c in cell:
    ATACFil = c + '_100bp.bedgraph'
    ATACSource = os.path.join(ATACFolder , ATACFil)
    ATACData = np.loadtxt(ATACSource , dtype=sig_type_1)
    ATACData_1 = np.loadtxt(os.path.join(ATACFolder , c + '_R1_100bp.bedgraph') , dtype=sig_type_1)
    ATACData_2 = np.loadtxt(os.path.join(ATACFolder , c + '_R2_100bp.bedgraph') , dtype=sig_type_1)
    for g in chrom:
        atac = ATACData[ATACData['chr'] == g]
        atac_1 = ATACData_1[ATACData_1['chr'] == g]
        atac_2 = ATACData_2[ATACData_2['chr'] == g]
        length = mm[g] // R
        new = np.zeros(length)
        new_1 = np.zeros(length)
        new_2 = np.zeros(length)
        for i in atac:
            start = i['start'] // R
            end = i['end'] // R
            score = i['score']
            n = end - start
            new[start] += score
            for j in np.arange(1,n):
                new[start + j] += score
        for i in new:
            all_data[c].append(i)
        for i in atac_1:
            start = i['start'] // R
            end = i['end'] // R
            score = i['score']
            n = end - start
            new_1[start] += score
            for j in np.arange(n):
                new_1[start + j] += score
        for i in new_1:
            all_data[c + '1'].append(i)
                
        for i in atac_2:
            start = i['start'] // R
            end = i['end'] // R
            score = i['score']
            n = end - start
            new_2[start] += score
            for j in np.arange(n):
                new_2[start + j] += score
                        
        for i in new_2:
            all_data[c + '2'].append(i)
                     

all_data_log2 = {'CCS' : [] , 'CCS1' : [] , 'CCS2':[] , 'NT5' : [] , 'NT51':[] , 'NT52':[] , 
            'NT6' : [] , 'NT61' : [] , 'NT62':[] , 'fESC' : [] , 'fESC1' : [] , 'fESC2' : []}

for c in all_data:
    for i in all_data[c]:
            all_data_log2[c].append(np.log2(i + 1))


all_cells = ['CCS_NT5' , 'CCS_NT6' , 'CCS_fESC' , 'NT5_NT6' , 'NT5_fESC' , 'NT6_fESC' , 
             'CCS1_CCS2' , 'NT51_NT52' , 'NT61_NT62' , 'fESC1_fESC2']


cor_data = {}
for c in all_cells:
    pp = PdfPages('D:\\Workspace_New\\Plot\\Correlation-Heatmap\\replicated_scatter\\ATAC\\'  + c + '_ATAC_log2_5K.pdf')
    c1 = c.split("_")[0]
    c2 = c.split("_")[1]
    cor = round(np.corrcoef(all_data_log2[c1] , all_data_log2[c2])[0][1],5)
    cor_data[c] = cor
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(all_data_log2[c1] , all_data_log2[c2] , alpha = 0.6 , c = 'red')
    if '1' in c:
        ax.set_xlabel(c1.split("1")[0] + '_R1' , size = 50)
        ax.set_ylabel(c2.split("2")[0] + '_R2' , size = 50)
    else:
        ax.set_xlabel(c1 , size = 50)
        ax.set_ylabel(c2 , size = 50)
    ax.set_xlim(0,14)
    ax.set_ylim(0,14)
    ax.text(2 , 12 , 'R = ' + str(cor) , size = 40 )
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
    pp.savefig(fig)
    pp.close()       

cor_matrix = np.zeros((4,4))

for i in cell:
    for j in cell:
        for k in all_cells:
            if (i in k) and (j in k):
                cor_matrix[cells[i]][cells[j]] = cor_data[k]
         
        
o = open('D:\\Workspace_New\\data\\correlation\\ATAC_correlation_matrix_log2(RPKM).txt' ,  'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()

            
            