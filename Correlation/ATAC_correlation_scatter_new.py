# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 15:46:39 2020

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
# import cPickle


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

cell = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
chrom=['1' , '2' , '3' , '4' , '5', '6' ,'7' ,'8' , '9' , '10' ,'11' , '12' ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X']
R = 100000

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['U4',np.int , np.int , np.float]})


ATACFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\signal\\normalization\\begraph_100bp'


f = open('E:\\Data\\literature_data\\genome\\mm10.txt' , 'r')
mm = {}
for i in f:
    i = i.strip().split()
    mm[i[0]] = int(i[1])
f.close()
    
#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6


all_data = {'CCS1' : [] , 'CCS2':[] , 'NT51':[] , 'NT52':[] , 'NT61' : [] , 'NT62':[] ,
            'F351':[] , 'F352':[] , 'F401' : [] , 'F402' : []}


for c in cell:
    ATACData_1 = np.loadtxt(os.path.join(ATACFolder , 'Uniq_' + c + '_R1.bedgraph') , dtype=sig_type_1)
    ATACData_2 = np.loadtxt(os.path.join(ATACFolder , 'Uniq_' + c + '_R2.bedgraph') , dtype=sig_type_1)
    for g in chrom:
        atac_1 = ATACData_1[ATACData_1['chr'] == g]
        atac_2 = ATACData_2[ATACData_2['chr'] == g]
        length = mm[g] // R + 1
        new_1 = np.zeros(length)
        new_2 = np.zeros(length)
        for i in atac_1:
            start = i['start'] // R
            score = i['score']
            new_1[start] += score
        all_data[c + '1'].extend(new_1)

        for i in atac_2:
            start = i['start'] // R
            score = i['score']
            new_2[start] += score
        all_data[c + '2'].extend(new_2)
                     

all_data_log2 = {}

for c in all_data:
    all_data_log2[c] = np.log2(np.array(all_data[c]) + 1)


all_cells = ['CCS1_CCS2' , 'NT51_NT52' , 'NT61_NT62' , 'F351_F352' , 'F401_F402']


cor_data = {}
for c in all_cells:
    pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\correlation\\'  + c + '_ATAC_log2_100K.pdf')
    c1 = c.split("_")[0]
    c2 = c.split("_")[1]
    cor = round(np.corrcoef(all_data_log2[c1] , all_data_log2[c2])[0][1],4)
    cor_data[c] = cor
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(all_data_log2[c1] , all_data_log2[c2] , alpha = 0.6 , c = 'red')
    ax.set_xlabel(c1.split("1")[0] + '_R1' , size = 50)
    ax.set_ylabel(c2.split("2")[0] + '_R2' , size = 50)
    
    ax.set_xlim(0,14)
    ax.set_ylim(0,14)
    ax.text(2 , 12 , 'R = ' + str(cor) , size = 40 )
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
    pp.savefig(fig)
    pp.close()       

      
 


           