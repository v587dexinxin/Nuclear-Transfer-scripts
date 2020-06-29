# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 18:58:08 2019

@author: xxli
"""


from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os

import matplotlib.pyplot as plt


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
                   labelleft = True, labelright = False , length = 5 ,labelsize = 30  )

def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
def getmatrix(inter,l_bin,r_bin):
    inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin),dtype = float )
    mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
           (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
    inter_extract = inter[mask]
    for i in inter_extract:
        if i['bin1'] != i['bin2']:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
        else:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
    return inter_matrix
    
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' ,'X']
cells = ['CCS','NT5','NT6','fESC']



CCS_Lib = np.load('D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\fESC_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT5_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT6_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}

CCS_tad = np.loadtxt('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\CCS_40K_allreps_union_small_domain.txt' , dtype = tad_type)
fESC_tad = np.loadtxt('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\fESC_40K_allreps_union_small_domain.txt' , dtype = tad_type)
NT5_tad = np.loadtxt('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\NT5_40K_allreps_union_small_domain.txt' , dtype = tad_type)
NT6_tad = np.loadtxt('D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\NT6_40K_allreps_union_small_domain.txt' , dtype = tad_type)
Tad = {'CCS':CCS_tad,
            'fESC':fESC_tad,
            'NT5':NT5_tad,
            'NT6':NT6_tad}

for c in cells:
    pp = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig3\\Boundary_nearby_interaction\\' + c + '_Boundary_nearby_interaction_40K.pdf')
    HiC_Lib = HiC_Data[c]
    Tad_Lib = Tad[c]
    matrix = np.zeros((40,40))
    for i in Tad_Lib:
        chro = i['chr']
        boundary = i['start']
        startHiC = (boundary - 800000) // 40000
        endHiC = (boundary + 800000) // 40000
        HiC_lib = HiC_Lib[chro]
        if endHiC > len(HiC_lib):
            continue
        matrix += HiC_lib[startHiC:endHiC , startHiC:endHiC]
    matrix = matrix / (len(Tad_Lib))
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 85)
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
    ticks = [0,20,40]
    labels = ['-800K','Boundary','+800K']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, rotation = 'horizontal')
    ax.set_xlabel(c , fontsize = 50 )
    caxis_H(ax)
    ##-----------------------colorbar--------------------------------------------------------------------------
    ax = fig.add_axes([Left + 0.5 , HB - 0.11 , 0.1 , 0.035])
    cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
    cbar.set_ticks([0 , int(vmax)])
    caxis_colorbar(ax)
    
    pp.savefig(fig)
    pp.close()
        