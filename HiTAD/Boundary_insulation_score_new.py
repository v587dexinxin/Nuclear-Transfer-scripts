# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:41:01 2020

@author: xxli
"""


from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import sys
import os

import matplotlib.pyplot as plt



 
def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new
                               

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

    
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S8' , np.int , np.int]})

boundary_type = np.dtype({'names':['chr' , 'pos'],
                     'formats':['S8' , np.int]})


chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
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




 

#-------------------------------------insulation_score-----------------------------------------           
Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/union_Filtered_Boundary_40K.txt' , dtype = boundary_type)
#Tad = Tad[Tad['chr'] != 'X']


for c in cells:
    out = open('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Insulation_score/' + c + '_Insulation_score_40K.txt' , 'w')
    out.writelines('\t'.join(['chr' , 'Boundary' , 'Insulation_Score']) + '\n')
    HiC_Lib = HiC_Data[c]
    Tad_Lib = Tad
    for i in Tad_Lib:
        chro = i['chr']
        boundary = i['pos']
        startHiC = (boundary - 600000) // 40000
        endHiC = (boundary + 600000) // 40000
        HiC_lib = HiC_Lib[chro]
        if endHiC > len(HiC_lib):
            continue
        matrix = HiC_lib[startHiC:endHiC , startHiC:endHiC]
        insulation = acquireSingleIns(matrix , 15 , 10 , 'divide')
        out.writelines('\t'.join([chro , str(boundary) , str(insulation)]) + '\n')
        
    out.close()
        

#-------------------------------------union_ins2onefil----------------------------------

boundary_type = np.dtype({'names' : ['chr' , 'site'] , 
                          'formats' : ['S8' , np.int]})
classify = ['CCS_noNT_nofESC' , 'CCS_noNT_fESC' , 'CCS_NT_nofESC' , 'noCCS_NT_fESC' , 
            'noCCS_NT_nofESC' , 'noCCS_noNT_fESC' , 'CCS_NT_fESC']

for cl in classify:              
    Tad = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/HiTAD/40K/respective/stable_600K/boundary_Venn3_classify/Boundary_' + cl + '.txt' , dtype = boundary_type)
    Tad = Tad[Tad['chr'] != 'X']
    out = open('/public/home/xxli/data/BDF1_New/HiC/HiTAD/40K/respective/stable_600K/boundary_Venn3_classify/Insulation_score/' + cl + '_Insulation_score_40K.txt' , 'w')
    out.writelines('\t'.join(['chr' , 'Boundary' , 'CCS_Ins' , 'NT5_Ins' , 'NT6_Ins' , 'fESC_Ins']) + '\n')
    
    for i in Tad:
        chro = i['chr']
        boundary = i['site']
        startHiC = (boundary - 300000) // 20000
        endHiC = (boundary + 300000) // 20000
        ins = []
        for c in cells:
            HiC_Lib = HiC_Data[c]
            Tad_Lib = Tad
            matrix = np.zeros((30,30))
            HiC_lib = HiC_Lib[chro]
            if endHiC > len(HiC_lib):
                continue
            matrix += HiC_lib[startHiC:endHiC , startHiC:endHiC]
            insulation = acquireSingleIns(matrix , 15 , 10 , 'divide')
            ins.append(insulation)
        
        out.writelines('\t'.join([chro , str(boundary) , str(ins[0]) , str(ins[1]) , str(ins[2]) , str(ins[3])]) + '\n')
        
    out.close()
    
    