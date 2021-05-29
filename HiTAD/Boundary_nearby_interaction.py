# -*- coding: utf-8 -*-
"""
Created on Sat Mar 09 10:19:01 2019

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
                                           ['#FFFFFF' , '#CD0000'])
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
        
    


def acquireSingleIns(matrix_data_chr,bound,left_right,category):
    ins=0
    start_site=0;end_site=matrix_data_chr.shape[0]
    if ((bound-left_right<start_site)|(end_site<bound+left_right)):        
        return ins    

    aa=matrix_data_chr[bound-left_right+1:bound+1,bound+1:bound+left_right+1]
    b1=[[matrix_data_chr[i,j] for i in range(bound-left_right,bound) if j>i] 
            for j in range(bound-left_right,bound)]
    b2=[[matrix_data_chr[i,j] for i in range(bound+1,bound+left_right+1) if j>i] 
            for j in range(bound+1,bound+left_right+1)]
    
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


def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S8' , np.int , np.int]})
boundary_type = np.dtype({'names':['chr' , 'pos'],
                     'formats':['S8' , np.int]})




cells = ['CCS','NT5','NT6','F35','F40']



#HiC Data Process
CCS_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/CCS_matrix/Cor_CCS_40K.npz' , allow_pickle=True)
CCS_Lib = Get_nan_zero_Matrix(CCS_Lib)
NT5_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT5_matrix/Cor_NT5_40K.npz' , allow_pickle=True)
NT5_Lib = Get_nan_zero_Matrix(NT5_Lib)
NT6_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT6_matrix/Cor_NT6_40K.npz' , allow_pickle=True)
NT6_Lib = Get_nan_zero_Matrix(NT6_Lib)
F35_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F35_matrix/Cor_F35_40K.npz' , allow_pickle=True)
F35_Lib = Get_nan_zero_Matrix(F35_Lib)
F40_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F40_matrix/Cor_F40_40K.npz' , allow_pickle=True)
F40_Lib = Get_nan_zero_Matrix(F40_Lib)

HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'F35':F35_Lib,
            'F40':F40_Lib}


Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/union_Filtered_Boundary_40K.txt' , dtype = boundary_type)
Tad = Tad[Tad['chr'] != 'X']

Matrix = {}
for c in cells:
    pp = PdfPages(c + '_Boundary_nearby_interaction_40K_ave.pdf')
    HiC_Lib = HiC_Data[c]
    Tad_Lib = Tad
    matrix = np.zeros((30,30))
    for i in Tad_Lib:
        chro = i['chr']
        boundary = i['pos']
        startHiC = (boundary - 600000) // 40000
        endHiC = (boundary + 600000) // 40000
        HiC_lib = HiC_Lib[chro]
        if endHiC > len(HiC_lib):
            continue
        matrix += HiC_lib[startHiC:endHiC , startHiC:endHiC] 
    matrix = matrix / len(Tad_Lib)
    Matrix[c] = matrix
    insulation = acquireSingleIns(matrix , 15 , 10 , 'divide')
    print (insulation)
    oe_matrix = OE_matrix(matrix)
    nonzero = oe_matrix[np.nonzero(oe_matrix)]
    vmin = nonzero.min()
    vmax = np.percentile(nonzero, 95)
    print (insulation)
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(oe_matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(oe_matrix), 0, len(oe_matrix)), vmax = 0.75, vmin = -0.75 , origin = 'lower')
    ticks = [0,15,30]
    labels = ['-600K','Boundary','+600K']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, rotation = 'horizontal')
    ax.set_xlabel(c + '_Insulation' + str(np.round(insulation,4)) + ', Tad_number:' + str(len(Tad_Lib)), fontsize = 20 )
    caxis_H(ax)
    ##-----------------------colorbar--------------------------------------------------------------------------
    ax = fig.add_axes([Left + 0.5 , HB - 0.11 , 0.1 , 0.035])
    cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
    cbar.set_ticks([np.round(vmin , 2) + 0.01 , np.round(vmax , 2) - 0.01])
    cbar.set_ticks([-0.75 , 0.75])
    caxis_colorbar(ax)
    
    pp.savefig(fig)
    pp.close()
        


