# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 18:37:45 2019

@author: xxli
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os


# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
                

def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 ,labelsize = 20  )

def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
def OE_matrix(matrix):
    matrix_new = np.zeros((len(matrix) , len(matrix)))
    exp = []
    for i in range(len(matrix)):
        tmp = []
        for j in range(len(matrix) - i):
            tmp.append((matrix[j][j + i]))
        m = sum(tmp)/len(tmp)
        if m == 0:
            exp.append(1)
        else:
            exp.append(m)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            dis = abs(i-j)
            if matrix[i][j] == 0:
                matrix_new[i][j] = 0
            else:
                matrix_new[i][j] = np.log2(matrix[i][j] / exp[dis])
    return matrix_new
        
def Get_Loop_Strength(matrix,start,end):
    """
    """
    
    matrix1 = matrix[start - 1 : start + 2 , end - 1 : end + 2]
    matrix2 = matrix[start + 2 : start + 5 , end + 2 : end + 5]
    
    n1 = len(np.nonzero(matrix1)[0])
    n2 = len(np.nonzero(matrix2)[0])
    ave1 = matrix1.sum() / n1
    ave2 = matrix2.sum() / n2
    strength = ave1 / ave2
    return strength    

                
size = (12, 12)
Left = 0.3 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

TadFolder = '/public/home/xxli/data/BDF1_New/HiC/HiTAD/respective/stable_600K'
cells = ['CCS','NT5','NT6','fESC']


CCS_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}

tads = np.loadtxt(os.path.join(TadFolder , 'union_domain.txt') , dtype = tad_type)

for c in cells:
    pp = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/average_interaction/' + c + '_tads_average_interaction_1.pdf')
    matrix = np.zeros((30, 30))
    for i in tads:
        g = i['chr']
        start = i['start']
        end = i['end']
        if i['end'] - i['start'] < 300000:
            continue
        center = (i['start'] + i['end']) / 2
        s = int((center - 300000) // 20000)
        e = int((center + 300000) // 20000)
        lib = HiC_Data[c][g]
        matrix += lib[s:e , s:e]
    matrix = matrix / (len(tads))
    loop_strength = Get_Loop_Strength(matrix , 8 , 23)
    print loop_strength
    oe_matrix = OE_matrix(matrix)
#    
#    for i in range(50):
#        for j in range(50):
#            if i == j:
#                matrix[i][j] = 0
#            else:
#                pass
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 80)
    print oe_matrix.min() , oe_matrix.max()
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(oe_matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(oe_matrix), 0, len(oe_matrix)), vmax = 0.5, vmin = -0.5 , origin = 'lower')
    ticks = [0 , 7.5 , 22.5 , 30]
    labels = ['-150K','Boundary','Boundary','+150K']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, rotation = 'horizontal')
    ax.set_xlabel(c + '_loop_' + str(np.round(loop_strength , 4 )), fontsize = 30 )
    caxis_H(ax)
    ##-----------------------colorbar--------------------------------------------------------------------------
    ax = fig.add_axes([Left + 0.5 , HB - 0.11 , 0.1 , 0.035])
    cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
    cbar.set_ticks([-0.5 , 0.5])
    caxis_colorbar(ax)
    
    pp.savefig(fig)
    pp.close()
    
    
    