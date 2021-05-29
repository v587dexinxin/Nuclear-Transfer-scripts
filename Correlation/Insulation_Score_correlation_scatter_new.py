# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 17:17:36 2021

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

    

def Tads2Boundary(Tads):
    boundary =[]
    for i in Tads:
        chro = i['chr']
        start = i['start']
        end = i['end']
        boundary.append((chro , start))
        boundary.append((chro , end))
    boundary = list(set(boundary))
    boundary.sort(key = lambda x:(x[0] , x[1]))
    boundary = np.array(boundary , dtype = boun_type)
    return boundary
    
    
    
def union_boundary(data_list):
   boundary_all = []
   union = []
   for c in data_list:
       for i in c:
           boundary_all.append(i)          
   boundary_all.sort(key = lambda x:(x[0] , x[1]))
   
   boun_number = len(boundary_all)
   for i in range(boun_number-1 , -1 , -1):
       site = boundary_all[i]['site']
       start = boundary_all[i-1]['site'] - 40000
       end = boundary_all[i-1]['site'] + 40000
       if (site >= start) and (site <= end):
           boundary_all.remove(boundary_all[i])
       else:
           pass

   union = np.array(boundary_all , dtype =  boun_type)
   

   return union
    
       
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6


tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S8' , np.int , np.int]})

boun_type = np.dtype({'names':['chr' , 'site' ],
                     'formats':['S8' , np.int ]})




chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
cells = ['CCS','NT5','NT6','F35','F40']



#HiC Data Process
CCS_R1_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/CCS_matrix/Cor_CCS_R1_40K.npz' , allow_pickle=True)
CCS_R1_Lib = Get_nan_zero_Matrix(CCS_R1_Lib)
CCS_R2_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/CCS_matrix/Cor_CCS_R2_40K.npz' , allow_pickle=True)
CCS_R2_Lib = Get_nan_zero_Matrix(CCS_R2_Lib)

NT5_R1_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT5_matrix/Cor_NT5_R1_40K.npz' , allow_pickle=True)
NT5_R1_Lib = Get_nan_zero_Matrix(NT5_R1_Lib)
NT5_R2_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT5_matrix/Cor_NT5_R2_40K.npz' , allow_pickle=True)
NT5_R2_Lib = Get_nan_zero_Matrix(NT5_R2_Lib)

NT6_R1_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT6_matrix/Cor_NT6_R1_40K.npz' , allow_pickle=True)
NT6_R1_Lib = Get_nan_zero_Matrix(NT6_R1_Lib)
NT6_R2_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT6_matrix/Cor_NT6_R2_40K.npz' , allow_pickle=True)
NT6_R2_Lib = Get_nan_zero_Matrix(NT6_R2_Lib)

F35_R1_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F35_matrix/Cor_F35_R1_40K.npz' , allow_pickle=True)
F35_R1_Lib = Get_nan_zero_Matrix(F35_R1_Lib)
F35_R2_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F35_matrix/Cor_F35_R2_40K.npz' , allow_pickle=True)
F35_R2_Lib = Get_nan_zero_Matrix(F35_R2_Lib)

F40_R1_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F40_matrix/Cor_F40_R1_40K.npz' , allow_pickle=True)
F40_R1_Lib = Get_nan_zero_Matrix(F40_R1_Lib)
F40_R2_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F40_matrix/Cor_F40_R2_40K.npz' , allow_pickle=True)
F40_R2_Lib = Get_nan_zero_Matrix(F40_R2_Lib)

HiC_Data = {'CCS': [CCS_R1_Lib, CCS_R2_Lib] , 
            'NT5': [NT5_R1_Lib, NT5_R2_Lib] , 
            'NT6': [NT6_R1_Lib, NT6_R2_Lib] , 
            'F35': [F35_R1_Lib, F35_R2_Lib] , 
            'F40': [F40_R1_Lib, F40_R2_Lib]}


#-------------------------------------insulation_score-----------------------------------------           

#Tad = Tad[Tad['chr'] != 'X']


pp = PdfPages('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Replicates/correlation_plot/Insulation_score_Correlation_replicates_scatter_40.pdf')

for c in cells:
    insulation_1 = [] ; insulation_2 = []
    Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Replicates/' + c + '_R1_R2_common_Domain_bottom_40K_respective_stable_600K.txt' , dtype = tad_type)
    Tad = Tad[Tad['chr'] != 'X']
    Tad_b_0 = Tads2Boundary(Tad)
    Tad_b = union_boundary([Tad_b_0])
    for i in Tad_b:
        chro = i['chr']
        boundary = i['site']
        startHiC = (boundary - 600000) // 40000
        endHiC = (boundary + 600000) // 40000
        HiC_lib_1 = HiC_Data[c][0][chro]
        HiC_lib_2 = HiC_Data[c][1][chro]
        if endHiC > len(HiC_lib_1):
            continue
        matrix1 = HiC_lib_1[startHiC:endHiC , startHiC:endHiC]
        matrix2 = HiC_lib_2[startHiC:endHiC , startHiC:endHiC]
        insulation1 = acquireSingleIns(matrix1 , 15 , 10 , 'divide')
        insulation2 = acquireSingleIns(matrix2 , 15 , 10 , 'divide')
        insulation_1.append(insulation1)
        insulation_2.append(insulation2)
    
    vmax1 = np.percentile(insulation_1 , 99.999)
    vmin1 = np.percentile(insulation_1 , 0.001)
    vmax2 = np.percentile(insulation_2 , 99.999)
    vmin2 = np.percentile(insulation_2 , 0.001)
    vmax = max(vmax1 , vmax2)
    vmin = min(vmin1 , vmin2)
    
    data1 = [] ; data2 = []
    for i in range(len(insulation_1)):
        if (insulation_1[i] > vmin) and (insulation_1[i] < vmax) and (insulation_2[i] > vmin) and (insulation_2[i] < vmax):
            data1.append(insulation_1[i])
            data2.append(insulation_2[i])
    
    cor = round(np.corrcoef(data1 , data2)[0][1],5)
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(data1 , data2 , alpha = 0.8 , c = 'red')
    # ax.set_xlim(-0.13,0.13)
    # ax.set_ylim(-0.13,0.13)
    # ax.set_xticks([-0.13 , 0.13])
    # ax.set_xticklabels([-0.13 , 0.13] , size = 25 )
    # ax.set_yticks([-0.13 , 0.13])
    # ax.set_yticklabels([-0.13 , 0.13] , size = 25 )
    ax.text(1.5 , 2.8 , 'PCC = ' + str(np.round(cor,3)) , size = 40 )
    ax.set_xlabel(c + '_R1' , size = 50 , labelpad = 0.001)
    ax.set_ylabel(c + '_R2' , size = 50 , labelpad = 0.001)
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
#    caxis_H(ax)
    pp.savefig(fig)
pp.close()
    
        
        
        
        
        
        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
