# -*- coding: utf-8 -*-
"""
Created on Tue Nov 05 12:17:41 2019

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


# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

Tad_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                    'formats':['S8' , np.int , np.int]})
                    
def normalize_TAD_interaction(Lib , Tad , length):           
    matrix = np.zeros((length , length))      
    n = 0   
    IS = np.zeros(length - 20)                                        
    for t in Tad:    
        chro = t['chr']
        start = t['start'] // R
        end = t['end'] // R
        if end - start < 20:
            continue
        point = {}
        for a in range(length * length):
            point[a] = []
    
        tad_length = end - start
        HiC_start = start - tad_length
        HiC_end = end + tad_length
        matrix0 = Lib[chro][HiC_start:HiC_end , HiC_start:HiC_end]
        if len(matrix0) < length:
            continue
        part = len(matrix0) / length
        for i in range(len(matrix0)):
            x = i // part
            for j in range(len(matrix0)):
                y = j // part
                site = x * length + y
                point[site].append(matrix0[i , j])
        tmp = np.zeros((length , length))
        for a in range(length * length):
            value = np.array(point[a]).mean()
            x1 = a // length
            y1 = a % length
            tmp[x1 , y1] += value
            
        for bound in range(10 , length - 10):
            Is = acquireSingleIns(tmp , bound , 10 , 'divide')
            IS[bound - 10] += Is
        matrix += tmp
        n += 1
    matrix = matrix / n        
    IS = IS / n
    print n
    return matrix , IS

        

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
    vmax = 0.001480
    vmin = 0.000022
    
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

    
#=============HeatMap + colorbar=================================
matrix_CCS , IS_CCS = normalize_TAD_interaction(CCS , Tad , 60)
matrix_NT5 , IS_NT5 = normalize_TAD_interaction(NT5 , Tad , 60)
matrix_NT6 , IS_NT6 = normalize_TAD_interaction(NT6 , Tad , 60)
matrix_fESC , IS_fESC = normalize_TAD_interaction(fESC , Tad , 60)

fig_CCS = Heatmap_Plot(matrix_CCS[10:50, 10:50] , 'CCS')
fig_NT5 = Heatmap_Plot(matrix_NT5[10:50, 10:50] , 'NT5')
fig_NT6 = Heatmap_Plot(matrix_NT6[10:50, 10:50] , 'NT6')
fig_fESC = Heatmap_Plot(matrix_fESC[10:50, 10:50] , 'fESC')

run_Plot(fig_CCS , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\CCS_average_Tad_interaction_40K.pdf')
run_Plot(fig_NT5 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\NT5_average_Tad_interaction_40K.pdf')
run_Plot(fig_NT6 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\NT6_average_Tad_interaction_40K.pdf')
run_Plot(fig_fESC , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\fESC_average_Tad_interaction_40K.pdf')


data = [IS_CCS , IS_NT5 , IS_NT6 , IS_fESC]
fig = Insulation_Plot(data)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\Normalized_Tad_interaction\\Average_Tad_interaction_40K.pdf')

