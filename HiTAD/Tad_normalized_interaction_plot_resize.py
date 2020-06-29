# -*- coding: utf-8 -*-
"""
Created on Wed Jun 03 19:33:47 2020

@author: han-luo
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
from skimage.transform import resize


# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

Tad_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                    'formats':['S8' , np.int , np.int]})

def normalize_TAD_interaction(Lib , Tad , length = 40):
    resize_matrix = np.zeros((length , length))
    n = 0
    for i in Tad:
        chro = i['chr']
        start = i['start'] // R 
        end = i['end'] // R
        d = end - start 
        if end - start < 40:
            continue
        startHiC = start - d // 2
        endHiC = end + d // 2
        matrix = Lib[chro][startHiC:endHiC , startHiC:endHiC]
        resize_m = resize(matrix , (length , length), order = 3)
        resize_matrix += resize_m 
        n += 1
    resize_matrix = resize_matrix / n
    return resize_matrix 
    
    
    
def acquireSingleIns(matrix_data_chr,bound,left_right,category):
    ins=0
    start_site=0;end_site=matrix_data_chr.shape[0]
    if ((bound-left_right<start_site)|(end_site<bound+left_right)):        
        return ins    
    aa=matrix_data_chr[bound-left_right:bound-0,bound+0:bound+left_right]
    b1=[[matrix_data_chr[i,j] for i in range(bound-left_right,bound-0) if j>i] 
            for j in range(bound-left_right,bound-0)]
    b2=[[matrix_data_chr[i,j] for i in range(bound+0,bound+left_right) if j>i] 
            for j in range(bound+0,bound+left_right)]
    
    aa_zero=sum([sum(np.array(item)==0) for item in aa])
    b1_zero=sum([sum(np.array(item)==0) for item in b1])
    b2_zero=sum([sum(np.array(item)==0) for item in b2])
    if aa_zero+b1_zero+b2_zero>=left_right:
        return ins    
    aa_sum=sum([sum(item) for item in aa])
    b1_sum=sum([sum(item) for item in b1])
    b2_sum=sum([sum(item) for item in b2])
    if aa_sum>0:
        if(category=='divide'):
            ins=np.log2((aa_sum+b1_sum+b2_sum)/float(aa_sum))
        elif(category=='average'):
            ins=aa_sum/float(left_right)/left_right
        else:
            print('the calc type went wrong')
    return ins
    
def Heatmap_Plot(matrix ,  c):
    matrix = matrix / matrix.sum()
    vmax = np.percentile(matrix , 90)
    vmin = np.percentile(matrix , 5)
    
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                    extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, vmin = vmin , origin = 'lower')
    cxlim = ax1.get_xlim()
    cylim = ax1.get_ylim()
    ## Ticks and Labels
    ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
    labels = ['-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD']
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_xlabel(c , fontsize = 30 )
    
    ax1.set_xlim(cxlim)
    ax1.set_ylim(cylim)                    
    ## Colorbar
    ax2 = fig.add_axes([Left + 0.6 , HB - 0.06 , 0.1 , 0.035])
    fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [vmin ,vmax])
    return fig

def Insulation_Plot(data):
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    ax1.plot(np.arange(len(data[0])) , data[0] , label = 'CCS')
    ax1.plot(np.arange(len(data[1])) , data[1] , label = 'NT5')
    ax1.plot(np.arange(len(data[2])) , data[2] , label = 'NT6')
    ax1.plot(np.arange(len(data[3])) , data[3] , label = 'fESC')
    
    
    ## Ticks and Labels
#    ticks = list(np.linspace(0 , len(matrix) , 7).astype(float))
#    labels = ['-TAD' , '-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD' , '+TAD']
#    ax1.set_xticks(ticks)
#    ax1.set_xticklabels(labels)
#    ax1.set_yticks(ticks)
#    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_ylabel('Insulation Score' , fontsize = 20 )
    ax1.legend()
    

    return fig

    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
        
    
R = 40000
size = (11, 12)
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7

##Load HiC_data
CCS = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
CCS = CCS['Matrix'][()]
NT5 = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT5_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT5 = NT5['Matrix'][()]
NT6 = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT6_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT6 = NT6['Matrix'][()]
fESC = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\fESC_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
fESC = fESC['Matrix'][()]

##Load Tad_data
Tad = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\union_domain.txt' , dtype = Tad_type)
Tad = Tad[Tad['chr'] != 'X']

matrix_CCS = normalize_TAD_interaction(CCS , Tad , 40)
matrix_NT5 = normalize_TAD_interaction(NT5 , Tad , 40)
matrix_NT6 = normalize_TAD_interaction(NT6 , Tad , 40)
matrix_fESC = normalize_TAD_interaction(fESC , Tad , 40)
Heatmap_Plot(matrix_CCS , 'CCS')
Heatmap_Plot(matrix_NT5 , 'NT5')
Heatmap_Plot(matrix_NT6 , 'NT6')
Heatmap_Plot(matrix_fESC , 'fESC')
                    