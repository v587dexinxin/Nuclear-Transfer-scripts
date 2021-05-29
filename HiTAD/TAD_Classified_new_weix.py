# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 19:05:47 2021

@author: xxli
"""

from __future__ import division
import numpy as np 
import pyBigWig
import pandas as pd
import scipy
from scipy import stats
import matplotlib 
# Use a non-interactive backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import wilcoxon
import math
import os
    
    
def Sig_To_100bp(signal , control):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        if control == 'input':
            tmp_data = np.array(list(signal.intervals(g)) , dtype = signal_type)
        else:
            tmp_data = np.array(list(signal.intervals(g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            # end = line['end'] // 100
            New_Data[g][start] += line['value']
    
    return New_Data
    

def Get_tads(TadFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    tad_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    tads = []
    fileSource = os.path.join(TadFil)
    Tad = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = tad_type)
    Tad = Tad[Tad['chr'] != 'X']
    
    for i in Tad:
        if i['end'] - i['start'] >= 400000:
            tads.append((i['chr'].lstrip('chr') , i['start'] , i['end']))
        else:
            continue
    tads = np.array(tads , dtype = tad_type)
    return tads
    
def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new



def calc_consolidation(a, b, contacts,cutoff): #contacts: {(start, end):value,...}
    if a > b: a,b = b,a
    D=b-a+1
    intraTAD,outTAD = [0.0,0.0],[0.0,0.0]
    for x in range(a-int(D/2),b+1):
        for y in range(a,b+int(D/2)+2):
            if x+cutoff <= y and y <= x+b-a and y+x <= 2*b and y+x >= 2*a:
                if x>=a and y<=b:
                    intraTAD[1] +=1
                    try:
                        intraTAD[0] += contacts[(x,y)]
                    except:
                        pass
                else:
                    outTAD[1] +=1
                    try:
                        outTAD[0] += contacts[(x,y)]
                    except:
                        pass
    
    try:
        score = (intraTAD[0]/intraTAD[1])/(outTAD[0]/outTAD[1])
    except:
        score = 0
    return intraTAD[0],intraTAD[1], outTAD[0], outTAD[1],score


def Box_plot(data , cl):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'orange','linewidth':1},
            medianprops={'color':'orange','linewidth':1},
            capprops={'color':'orange','linewidth':1},
            whiskerprops={'color':'orange','linewidth':1})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'green','linewidth':1},
            medianprops={'color':'green','linewidth':1},
            capprops={'color':'green','linewidth':1},
            whiskerprops={'color':'green','linewidth':1})
    ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':1},
            medianprops={'color':'fuchsia','linewidth':1},
            capprops={'color':'fuchsia','linewidth':1},
            whiskerprops={'color':'fuchsia','linewidth':1})
    
    ax.plot([0.5,5.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 , 4 , 5])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F35' , 'F40' ] , fontsize = 28)
    ax.set_xlim((0.5 , 5.5))
    ax.set_title(cl + '_numbers' + str(len(data[0])))
    return fig





def Box_plot_3cellline(data , cl):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'orange','linewidth':1},
            medianprops={'color':'orange','linewidth':1},
            capprops={'color':'orange','linewidth':1},
            whiskerprops={'color':'orange','linewidth':1})

    d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels(['CCs' , 'NTs' , 'fESC'] , fontsize = 10)
    ax.set_xlabel('CCS_NTs: ' + str(d1) + ' , CCS_fESC: ' + str(d2) + ' , NTs_fESC: ' + str(d3))
    ax.set_xlim((0.5 , 3.5))
    ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((-0.25 , 6))
    
    return fig

def Box_plot_3cellline_merged(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot([data[0] , data[3] , data[6] , data[9]] , positions=[1 , 5 , 9, 13] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot([data[1] , data[4] , data[7] , data[10]] , positions=[2 , 6 , 10, 14] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})
    ax.boxplot([data[2] , data[5] , data[8] , data[11]] , positions=[3 , 7 , 11, 15] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'orange','linewidth':1},
            medianprops={'color':'orange','linewidth':1},
            capprops={'color':'orange','linewidth':1},
            whiskerprops={'color':'orange','linewidth':1})

    d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    d4 = np.round(wilcoxon(data[3] , data[4])[1] , 5)
    d5 = np.round(wilcoxon(data[3] , data[5])[1] , 5)
    d6 = np.round(wilcoxon(data[4] , data[5])[1] , 5)
    
    d7 = np.round(wilcoxon(data[6] , data[7])[1] , 5)
    d8 = np.round(wilcoxon(data[6] , data[8])[1] , 5)
    d9 = np.round(wilcoxon(data[7] , data[8])[1] , 5)
    
    d10 = np.round(wilcoxon(data[9] , data[10])[1] , 5)
    d11 = np.round(wilcoxon(data[9] , data[11])[1] , 5)
    d12 = np.round(wilcoxon(data[10] , data[11])[1] , 5)

    
    ax.plot([4 , 4] ,  [-0.25 , 9] , linewidth = 1.5 , linestyle = ':' , color = 'grey')
    ax.plot([8 , 8] ,  [-0.25 , 9] , linewidth = 1.5 , linestyle = ':' , color = 'grey')
    ax.plot([12 , 12] ,  [-0.25 , 9] , linewidth = 1.5 , linestyle = ':' , color = 'grey')
    ax.set_xticks([1 , 2 , 3 , 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15, 16])
    ax.set_xticklabels(['CCs' , 'NTs' , 'fESC' , 'Repro' , 'CCs' , 'NTs' , 'fESC' , 'Resis' , 'CCs' , 'NTs' , 'fESC' , 'Over' , 'CCs' , 'NTs' , 'fESC' , 'Others'] , fontsize = 10)
    ax.set_xlabel('CCS_NTs: ' + str(d1) + ' , CCS_fESC: ' + str(d2) + ' , NTs_fESC: ' + str(d3) + '\n' + \
                  'CCS_NTs: ' + str(d4) + ' , CCS_fESC: ' + str(d5) + ' , NTs_fESC: ' + str(d6) + '\n' + \
                  'CCS_NTs: ' + str(d7) + ' , CCS_fESC: ' + str(d8) + ' , NTs_fESC: ' + str(d9) + '\n' + \
                  'CCS_NTs: ' + str(d10) + ' , CCS_fESC: ' + str(d11) + ' , NTs_fESC: ' + str(d12))
    ax.set_xlim((0.5 , 16.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((-0.25 , 6))
    
    return fig


def Box_plot_2cellline(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})


    d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)


    
    ax.set_xticks([1 , 2])
    ax.set_xticklabels(['Repro' , 'unRepro'] , fontsize = 10)
    ax.set_xlabel('Repro_unRepro: ' + str(d1))
    ax.set_xlim((0.5 , 2.5))
    # ax.set_ylim((0 , 25))
    
    return fig

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




def distance(sites_1,sites_2):
    dis = math.sqrt((sites_1[1] / 40000 - sites_2[1] / 40000)**2 + (sites_1[2] / 40000-sites_2[2] / 40000)**2)
    return dis


def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a



def Get_diff_loops(loop_1 , loop_2):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        for j in tmp:
            d = distance([chro , i[1] , i[2]] , j)
            if d <= initial_distance:
                common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        for j in tmp:
            d = distance([chro , i[1] , i[2]] , j)
            if d <= initial_distance:
                common_2.append(i)
    common_1 = np.array(common_1 , dtype = loop_type)
    common_2 = np.array(common_2 , dtype = loop_type)
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print (len(loop_1) , len(loop_2) , len(common_1) , len(common_2))
    for i in common_1:
        try:
            loop_1.remove(i)
        except:
            pass
    for j in common_2:
        try:
            loop_2.remove(j)
        except:
            pass
    print (len(loop_1) , len(loop_2))   
    diff_1 = np.array(Sort(loop_1 , 0 , 1 ), dtype = loop_type)
    diff_2 = np.array(Sort(loop_2 , 0 , 1 ), dtype = loop_type)
    return diff_1 , diff_2 , common_1 , common_2


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    
    
    
    
bound_type = np.dtype({'names':['chr' , 'pos'] , 
                    'formats':['S8' , np.int]})
tad_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                       'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})
union_type = np.dtype({'names':['chr' , 'start' , 'end' , 'CCS' , 'NTs' , 'fESC'],
                       'formats':['S8' , np.int , np.int , np.float , np.float , np.float]})


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
cutoff = 400000
resolution = 40000
R = 40000
initial_distance = math.sqrt(2) + 0.05
cells = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']

##-------------------------------tad_data--------------------------

union_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/union_domain_NTs_fESC_merged.txt')
union_tad = union_tad[union_tad['chr'] != 'X']
CCS_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/CCS_Domain_bottom_40K_respective_stable_400K.txt')
CCS_tad = CCS_tad[CCS_tad['chr'] != 'X']
NT5_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/NT5_Domain_bottom_40K_respective_stable_400K.txt')
NT5_tad = NT5_tad[NT5_tad['chr'] != 'X']
NT6_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/NT6_Domain_bottom_40K_respective_stable_400K.txt')
NT6_tad = NT6_tad[NT6_tad['chr'] != 'X']
F35_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/F35_Domain_bottom_40K_respective_stable_400K.txt')
F35_tad = F35_tad[F35_tad['chr'] != 'X']
F40_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/F40_Domain_bottom_40K_respective_stable_400K.txt')
F40_tad = F40_tad[F40_tad['chr'] != 'X']


NTs_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/NTs_Domain_bottom_40K_respective_stable_400K.txt')
NTs_tad = NTs_tad[NTs_tad['chr'] != 'X']
fESC_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/fESC_Domain_bottom_40K_respective_stable_400K.txt')
fESC_tad = fESC_tad[fESC_tad['chr'] != 'X']


# NT5_1 , NT6_1 , common_NT , common_NT_2 = Get_diff_loops(NT5_tad , NT6_tad)
# F35_1 , F40_1 , common_fESC , common_fESC_2 = Get_diff_loops(F35_tad , F40_tad)

# NTs_tad = common_NT
# fESC_tad = common_fESC

#----------------------------HiC Data Process--------------------------------------
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

fESC_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/fESC_matrix/Cor_fESC_40K.npz' , allow_pickle=True)
fESC_Lib = Get_nan_zero_Matrix(fESC_Lib)
NTs_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NTs_matrix/Cor_NTs_40K.npz' , allow_pickle=True)
NTs_Lib = Get_nan_zero_Matrix(NTs_Lib)


HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'F35':F35_Lib,
            'F40':F40_Lib,
            'NTs':NTs_Lib,
            'fESC':fESC_Lib}




##--------------------------Tads_classify_new--------------------------------



n = 0
union_RTI = []
for i in union_tad:
    g = i['chr']
    start = i['start'] // resolution
    end = i['end'] // resolution
    if end - start < 10:
        n += 1
        continue
    u , v , w , x , ccs_score = calc_consolidation(start, end , HiC_Data['CCS'][g] ,cutoff // resolution)
    u , v , w , x , nts_score = calc_consolidation(start, end , HiC_Data['NTs'][g] ,cutoff // resolution)
    u , v , w , x , fesc_score = calc_consolidation(start, end , HiC_Data['fESC'][g] ,cutoff // resolution)

    if np.isnan(ccs_score):
        ccs_score = 0
        # print i
    if np.isnan(nts_score):
        nts_score = 0
        # print i
    if np.isnan(fesc_score):
        fesc_score = 0
        # print i
    if np.isinf(ccs_score):
        continue
    if np.isinf(nts_score):
        continue
    if np.isinf(fesc_score):
        continue

    union_RTI.append((g , i['start'] , i['end'] , ccs_score , nts_score , fesc_score))


union_RTI = np.array(union_RTI , dtype = union_type)

##Specific
CCS_speci = []
for i in union_RTI:
    g = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = CCS_tad[CCS_tad['chr'] == g]
    tmp2 = fESC_tad[fESC_tad['chr'] == g]
    tmp3 = NTs_tad[NTs_tad['chr'] == g]
    mask1 = (tmp1['start'] >= (start - resolution)) & (tmp1['start'] <= (start + resolution)) & (tmp1['end'] >= (end - resolution)) & (tmp1['end'] <= (end + resolution))
    mask2 = (tmp2['start'] >= (start - resolution)) & (tmp2['start'] <= (start + resolution)) & (tmp2['end'] >= (end - resolution)) & (tmp2['end'] <= (end + resolution))
    mask3 = (tmp3['start'] >= (start - resolution)) & (tmp3['start'] <= (start + resolution)) & (tmp3['end'] >= (end - resolution)) & (tmp3['end'] <= (end + resolution))
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]
    overlap3 = tmp3[mask3]
    if (overlap1.size != 0) and (overlap2.size == 0) and ((overlap3.size == 0)) and (i['CCS'] - i['fESC'] > 0.2):
        CCS_speci.append(i)
CCS_speci = np.array(CCS_speci , dtype = union_RTI.dtype)

fESC_speci = []
for i in union_RTI:
    g = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = CCS_tad[CCS_tad['chr'] == g]
    tmp2 = fESC_tad[fESC_tad['chr'] == g]
    tmp3 = NTs_tad[NTs_tad['chr'] == g]
    mask1 = (tmp1['start'] >= (start - resolution)) & (tmp1['start'] <= (start + resolution)) & (tmp1['end'] >= (end - resolution)) & (tmp1['end'] <= (end + resolution))
    mask2 = (tmp2['start'] >= (start - resolution)) & (tmp2['start'] <= (start + resolution)) & (tmp2['end'] >= (end - resolution)) & (tmp2['end'] <= (end + resolution))
    mask3 = (tmp3['start'] >= (start - resolution)) & (tmp3['start'] <= (start + resolution)) & (tmp3['end'] >= (end - resolution)) & (tmp3['end'] <= (end + resolution))
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]
    overlap3 = tmp3[mask3]
    
    if (overlap1.size == 0) and (overlap2.size != 0) and ((overlap3.size != 0)) and (i['CCS'] - i['fESC'] < -0.2):
        fESC_speci.append(i)
fESC_speci = np.array(fESC_speci , dtype = union_RTI.dtype)



## Classify
A_B_Repro = [] ; A_B_Resis = [] ; A_B_Over = []
B_A_Repro = [] ; B_A_Resis = [] ; B_A_Over = [] ; Others = []


for i in union_RTI:
    g = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = CCS_tad[CCS_tad['chr'] == g]
    tmp2 = fESC_tad[fESC_tad['chr'] == g]
    mask1 = (tmp1['start'] >= (start - resolution)) & (tmp1['start'] <= (start + resolution)) & (tmp1['end'] >= (end - resolution)) & (tmp1['end'] <= (end + resolution))
    mask2 = (tmp2['start'] >= (start - resolution)) & (tmp2['start'] <= (start + resolution)) & (tmp2['end'] >= (end - resolution)) & (tmp2['end'] <= (end + resolution))
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]
    if ((overlap1.size != 0) and (overlap2.size == 0) and (i['CCS'] - i['fESC'] > 0.2 )):
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 1.2:
            A_B_Over.append(i)
        elif ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 0.7) and ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 1.2):
            A_B_Repro.append(i)
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 0.7:
            A_B_Resis.append(i)
    elif (overlap1.size == 0) and (overlap2.size != 0) and (i['CCS'] - i['fESC'] < -0.2 ):
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 1.2:
            B_A_Over.append(i)
        elif ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 0.7) and ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 1.2):
            B_A_Repro.append(i)
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 0.7:
            B_A_Resis.append(i)
    elif (i['CCS'] - i['fESC'] > -0.2) and (i['CCS'] - i['fESC'] < 0.2):
        Others.append(i)




A_B_Repro = np.array(A_B_Repro , dtype = union_RTI.dtype)
A_B_Resis = np.array(A_B_Resis , dtype = union_RTI.dtype)
A_B_Over = np.array(A_B_Over , dtype = union_RTI.dtype)

B_A_Repro = np.array(B_A_Repro , dtype = union_RTI.dtype)
B_A_Resis = np.array(B_A_Resis , dtype = union_RTI.dtype)
B_A_Over = np.array(B_A_Over , dtype = union_RTI.dtype)
Others = np.array(Others , dtype = union_RTI.dtype)


tads = {'A_B_Repro':A_B_Repro , 'A_B_Resis':A_B_Resis , 'A_B_Over':A_B_Over ,
        'B_A_Repro':B_A_Repro , 'B_A_Resis':B_A_Resis , 'B_A_Over':B_A_Over ,
        'Others':Others, 
        'Repro':np.hstack((A_B_Repro , B_A_Repro)) , 'Resis':np.hstack((A_B_Resis , B_A_Resis)) , 
        'Over':np.hstack((A_B_Over , B_A_Over)) , 
        'CCS_Speci':CCS_speci ,'fESC_Speci':fESC_speci}
     

##---------------------plot_3cellines--------------------------------    
    

for cl in tads :
    
    out = open('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/excel/' + cl + '_TADs.txt' , 'w')
    for i in tads[cl]:
        out.writelines('\t'.join([i['chr'] , str(i['start']) , str(i['end']) , str(i['CCS']) , str(i['NTs']) , str(i['fESC'])]) + '\n')
    out.close()
    # print cl , len(tads[cl])
    
    fig = Box_plot_3cellline([tads[cl]['CCS'] , tads[cl]['NTs'] , tads[cl]['fESC']] , cl)
    run_Plot(fig , '/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/plot/' + cl + '_Tad_Relative_TAD_intensity_3celllines_1.pdf')
    
    print cl , len(tads[cl])

    plt.close()
    

    
fig = Box_plot_3cellline_merged([tads['Repro']['CCS'] , tads['Repro']['NTs'] , tads['Repro']['fESC'] , tads['Resis']['CCS'] , tads['Resis']['NTs'] , tads['Resis']['fESC'] , tads['Over']['CCS'] , tads['Over']['NTs'] , tads['Over']['fESC'] , tads['Others']['CCS'] , tads['Others']['NTs'] , tads['Others']['fESC']])
run_Plot(fig , '/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_new/plot/Tad_Relative_TAD_intensity_3celllines_allclassify.pdf')
plt.close()


##---------------------heatmap_plot-------------------------------------


matrix = np.zeros((len(A_B_Repro) + len(B_A_Repro) + len(A_B_Resis) + len(B_A_Resis) + len(A_B_Over) + len(B_A_Over) + len(Others) , 3))














