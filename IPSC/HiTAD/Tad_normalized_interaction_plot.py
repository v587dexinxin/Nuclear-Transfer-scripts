# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:15:13 2019

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
    vmax = 0.0008
    vmin = 0
    
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
    ax1.plot(np.arange(len(data[0])) , data[0] , label = 'MEF')
    ax1.plot(np.arange(len(data[1])) , data[1] , label = 'IPS_P3')
    ax1.plot(np.arange(len(data[2])) , data[2] , label = 'E14')

    
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
        
chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']    
R = 40000
size = (11, 12)
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7

##Load HiC_data
MEF_0 = np.load('H:\\Workspace_New\\data\\IPSC\\HiC\\Matrix\\Corrected_MEF_matrix_40K.npz')
IPS_P3_0 = np.load('H:\\Workspace_New\\data\\IPSC\\HiC\\Matrix\\Corrected_IPS_P3_matrix_40K.npz')
E14_0 = np.load('H:\\Workspace_New\\data\\IPSC\\HiC\\Matrix\\Corrected_E14_matrix_40K.npz')

MEF = {} ; IPS_P3 = {} ; E14 = {}
for g in chroms:
#    MEF[g] = np.nan_to_num(MEF_0[g])
#    IPS_P3[g] = np.nan_to_num(IPS_P3_0[g])
    E14[g] = np.nan_to_num(E14_0[g])


##Load Tad_data
Tad = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\HiTAD\\MEF_TADs_40K\\MEF_TADs_40K_Domain_40K.txt' , dtype = Tad_type)
Tad = Tad[Tad['chr'] != 'X']

    
#=============HeatMap + colorbar=================================
matrix_MEF , IS_MEF = normalize_TAD_interaction(MEF , Tad , 60)
matrix_IPS_P3 , IS_IPS_P3 = normalize_TAD_interaction(IPS_P3 , Tad , 60)
matrix_E14 , IS_E14 = normalize_TAD_interaction(E14 , Tad , 60)

fig_MEF = Heatmap_Plot(matrix_MEF[10:50, 10:50] , 'MEF')
fig_IPS_P3 = Heatmap_Plot(matrix_IPS_P3[10:50, 10:50] , 'IPS_P3')
fig_E14 = Heatmap_Plot(matrix_E14[10:50, 10:50] , 'E14')


run_Plot(fig_MEF , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\MEF_average_Tad_interaction_40K.pdf')
run_Plot(fig_IPS_P3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPS_P3_average_Tad_interaction_40K.pdf')
run_Plot(fig_E14 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\E14_average_Tad_interaction_40K.pdf')


data = [IS_MEF , IS_IPS_P3 , IS_E14]
fig = Insulation_Plot(data)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\IPSC_Average_Tad_interaction_40K.pdf')

