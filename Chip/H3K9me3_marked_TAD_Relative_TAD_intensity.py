# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:36:54 2021

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
import os
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind


def Get_tads(TadFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    tad_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
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
    
    


 

def Sig_To_100bp(signal , chro):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals(chro + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            # end = line['end'] // 100
            # for i in range(start,end):
            New_Data[g][start] += line['value']
    
    return New_Data

    

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


def Box_plot(data , cl , percent):                
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
    ax.set_title(cl + '_numbers' + str(len(data[0])) + '_percent:' + str(percent))
    return fig


def Box_plot_3cellline(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkgreen','linewidth':1},
            medianprops={'color':'darkgreen','linewidth':1},
            capprops={'color':'darkgreen','linewidth':1},
            whiskerprops={'color':'darkgreen','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'mediumvioletred','linewidth':1},
            medianprops={'color':'mediumvioletred','linewidth':1},
            capprops={'color':'mediumvioletred','linewidth':1},
            whiskerprops={'color':'mediumvioletred','linewidth':1})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})

    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)
    d1 = np.round(ttest_rel(data[0] , data[1])[1] , 5)
    d2 = np.round(ttest_rel(data[0] , data[2])[1] , 5)
    d3 = np.round(ttest_rel(data[1] , data[2])[1] , 5)
    
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels(['CCs' , 'NTs' , 'fESC'] , fontsize = 10)
    ax.set_xlabel('CCS_NTs: ' + str(d1) + ' , CCS_fESC: ' + str(d2) + ' , NTs_fESC: ' + str(d3))
    ax.set_xlim((0.5 , 3.5))
    # ax.set_ylim((0 , 6.5))
    
    return fig


def Box_plot_4cellline(data):                
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
    
    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['Repro' , 'Resis' , 'Hyper' , 'Others'] , fontsize = 28)
    ax.set_xlim((0.5 , 4.5))
    # ax.set_title(cl + '_numbers' + str(len(data[0])) + '_percent:' + str(percent))
    return fig




def Box_plot_2cellline(data):                
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


    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    d1 = np.round(ttest_ind(data[0] , data[1])[1] , 5)

    
    ax.set_xticks([1 , 2])
    ax.set_xticklabels(['Repro' , 'abnormal'] , fontsize = 10)
    ax.set_xlabel('Repro_abnormal: ' + str(d1))
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

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
cutoff = 400000
res = 40000
cells = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']

##-------------------------------tad_data--------------------------
Repro_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/Repro_TADs.txt' , dtype = tad_type)
Resis_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/Resis_TADs.txt' , dtype = tad_type)
Over_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/Over_TADs.txt' , dtype = tad_type)
Others_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/Others_TADs.txt' , dtype = tad_type)
CCS_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/CCS_Speci_TADs.txt' , dtype = tad_type)
fESC_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/fESC_Speci_TADs.txt' , dtype = tad_type)
CCS_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/CCS_Domain_bottom_40K_respective_stable_400K.txt' , dtype = tad_type)
NTs_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/NTs_Domain_bottom_40K_respective_stable_400K.txt' , dtype = tad_type)
fESC_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/fESC_Domain_bottom_40K_respective_stable_400K.txt' , dtype = tad_type)
union_tad = Get_tads('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/union_domain_NTs_fESC_merged.txt')
union_tad = union_tad[union_tad['chr'] != 'X']

union_tad = np.hstack((Repro_Tad , Resis_Tad , Over_Tad , Others_Tad))

union_1 = []
for i in union_tad:
    chro = i['chr']
    start_s = i['start'] - 40000
    start_e = i['start'] + 40000
    end_s = i['end'] - 40000
    end_e = i['end'] + 40000  
    
    tmp1 = CCS_Speci_Tad[CCS_Speci_Tad['chr'] == chro]
    tmp2 = fESC_Speci_Tad[fESC_Speci_Tad['chr'] == chro]
    mask1 = (tmp1['start'] >=  start_s) & (tmp1['start'] <=  start_e) & (tmp1['end'] >=  end_s) & (tmp1['end'] <=  end_e)
    mask2 = (tmp2['start'] >=  start_s) & (tmp2['start'] <=  start_e) & (tmp2['end'] >=  end_s) & (tmp2['end'] <=  end_e)
    overlap1 = tmp1[mask1]
    overlap2 = tmp2[mask2]
    if (overlap1.size == 0) and (overlap2.size == 0):
        union_1.append(i)
union_1 = np.array(union_1 , dtype = union_tad.dtype)
    
    
Tad_data = {'Repro':Repro_Tad,
            'Resis':Resis_Tad,
            'Over':Over_Tad,
            'Others':Others_Tad,
            'CCS_Speci':CCS_Speci_Tad,
            'fESC_Speci':fESC_Speci_Tad,
            'CCS':CCS_Tad,
            'NTs':NTs_Tad,
            'fESC':fESC_Tad,
            'union_CCS':union_1,
            'union_fESC':union_1}


##-------------------------------chip_data----------------------------------------
# chip1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_chip_50bp.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_input_50bp.bw")
chip1 = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw')
input1 = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw')
chip2 = pyBigWig.open("/public/home/lixinxin/data/literature/cell_cooperative/mapping/Unique_ESC_H3K9me3_50bp.bw")
input2 = pyBigWig.open("/public/home/lixinxin/data/literature/cell_cooperative/mapping/Unique_ESC_H3K9me3_Input_50bp.bw")


Chip1 = Sig_To_100bp(chip1 , 'chr')
Input1 = Sig_To_100bp(input1 , 'chr')
Chip2 = Sig_To_100bp(chip2 , '')
Input2 = Sig_To_100bp(input2 , '')


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






    
# ##-----------------------Tads_classify_Relative TAD intensity (RTI)------------------------------------   

    

for cl in ['CCS_Speci' , 'fESC_Speci' , 'union_CCS' , 'union_fESC']:
    
    H3K9mes_tads = []
    
    if cl == 'CCS_Speci' or cl == 'union_CCS':
        out = open('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_TADs/' + cl + '_CCS_H3K9me3_marked_TADs_fc1.4.txt' , 'w')
        for g in chroms:
            tmp_tad = Tad_data[cl][Tad_data[cl]['chr'] == g]
            tmp1 = Chip1[g]
            tmp3 = Input1[g]
            sum1 = tmp1.sum()
            sum3 = tmp3.sum()
            for i in tmp_tad:
                start = i[1] // 100
                end = i[2] // 100
                data1 = tmp1[start:end].sum() / (end - start)
                data3 = tmp3[start:end].sum() / (end - start)
                if data1 > 1.4 * data3:
                    H3K9mes_tads.append(i)
    else:
        out = open('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_TADs/' + cl + '_fESC_H3K9me3_marked_TADs_fc1.4.txt' , 'w')
        for g in chroms:
            tmp_tad = Tad_data[cl][Tad_data[cl]['chr'] == g]
            tmp1 = Chip2[g]
            tmp3 = Input2[g]
            sum1 = tmp1.sum()
            sum3 = tmp3.sum()
            for i in tmp_tad:
                start = i[1] // 100
                end = i[2] // 100
                data1 = tmp1[start:end].sum() / (end - start)
                data3 = tmp3[start:end].sum() / (end - start)
                if data1 > 1.4 * data3:
                    H3K9mes_tads.append(i)        
          
            
    print cl , (len(H3K9mes_tads) , len(H3K9mes_tads) / len(Tad_data[cl]))
    for i in H3K9mes_tads:
        out.writelines('\t'.join([i[0] , str(i[1]) , str(i[2])]) + '\n')
    out.close()
    
    
    
    # Score = {}
    # for c in cells:
    #     Score[c] = []
    #     for i in H3K9mes_tads:
    #         g = i[0]
    #         start = i[1] // res
    #         end = i[2] // res
    #         if end - start <= 10:
    #             continue
    #         u , v , w , x , score = calc_consolidation(start, end , HiC_Data[c][g] ,cutoff // res)
    #         if np.isnan(score):
    #             score = 0
    #         Score[c].append(score)
    
 
    # fig = Box_plot([Score['CCS'] , Score['NT5'] , Score['NT6'] , Score['F35'] , Score['F40']] , cl , (len(H3K9mes_tads) / len(Tad_data[cl])))
    # run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/' + cl + '_H3K9me3_marked_Tad_Relative_TAD_intensity_2.pdf')

    # plt.close()
 
    
    Score = {}
    for c in ['CCS' , 'NTs' , 'fESC']:
        Score[c] = []
        for i in H3K9mes_tads:
            g = i[0]
            start = i[1] // res
            end = i[2] // res
            if end - start <= 10:
                continue
            u , v , w , x , score = calc_consolidation(start, end , HiC_Data[c][g] ,cutoff // res)
            if np.isnan(score):
                score = 0
            
            Score[c].append(score)
            
    Score_1 = {'CCS':[] , 'NTs':[] , 'fESC':[]}
    for i in range(len(Score['CCS'])):
        if np.isinf(Score['CCS'][i]) or np.isinf(Score['NTs'][i]) or np.isinf(Score['fESC'][i]):
            continue
        else:
            Score_1['CCS'].append(Score['CCS'][i])
            Score_1['NTs'].append(Score['NTs'][i])
            Score_1['fESC'].append(Score['fESC'][i])
            
    
    data = [Score_1['CCS'] , Score_1['NTs'] , Score_1['fESC']]
    # print np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # print np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # print np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)
    
    # print np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # print np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # print np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    print ttest_rel(data[0] , data[1])[1]
    print ttest_rel(data[0] , data[2])[1]
    print ttest_rel(data[1] , data[2])[1]
    

    # print ttest_ind(data[0] , data[1])[1]
    # print ttest_ind(data[0] , data[2])[1]
    # print ttest_ind(data[1] , data[2])[1]
    
    
    fig = Box_plot_3cellline([Score['CCS'] , Score['NTs'] , Score['fESC']])
    run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/' + cl + '_H3K9me3_marked_Tad_Relative_TAD_intensity_3cellines_7.pdf')
    
        

    

##-----------------------Tads_classify_Relative TAD intensity (RTI)_CCS_fESC------------------------------------   

    

for cl in ['CCS' , 'fESC']:
    H3K9mes_tads = []
    for g in chroms:
        tmp_tad = Tad_data[cl][Tad_data[cl]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_tad:
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() / (end - start)
            data3 = tmp3[start:end].sum() / (end - start)
            if data1 > data3:
                H3K9mes_tads.append(i)
          
            
    print cl , (len(H3K9mes_tads) / len(Tad_data[cl]))
    
    
    
    Score = {}
    for c in cells:
        Score[c] = []
        for i in H3K9mes_tads:
            g = i[0]
            start = i[1] // res
            end = i[2] // res
            if end - start <= 10:
                continue
            u , v , w , x , score = calc_consolidation(start, end , HiC_Data[c][g] ,cutoff // res)
            if np.isnan(score):
                score = 0
            Score[c].append(score)
    
 
    fig = Box_plot([Score['CCS'] , Score['NT5'] , Score['NT6'] , Score['F35'] , Score['F40']] , cl , (len(H3K9mes_tads) / len(Tad_data[cl])))
    run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/' + cl + '_H3K9me3_marked_Tad_Relative_TAD_intensity.pdf')

    plt.close()
 
    
    if cl in ['CCS' , 'fESC']:
        Score = {}
        for c in ['CCS' , 'NTs' , 'fESC']:
            Score[c] = []
            for i in H3K9mes_tads:
                g = i[0]
                start = i[1] // res
                end = i[2] // res
                if end - start <= 10:
                    continue
                u , v , w , x , score = calc_consolidation(start, end , HiC_Data[c][g] ,cutoff // res)
                if np.isnan(score):
                    score = 0
                Score[c].append(score)
        
        
        fig = Box_plot_3cellline([Score['CCS'] , Score['NTs'] , Score['fESC']])
        run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/' + cl + '_H3K9me3_marked_Tad_Relative_TAD_intensity_3cellines.pdf')
        
        

    




##-----------------------Tads_classify_Chip_signal_4classifies------------------------------------   

# tads = {'Repro':Tad_data['Repro'] , 'un_Repro':np.hstack((Tad_data['Resis'] , Tad_data['Over'] , Tad_data['Others']))}
tads = {'Repro':Tad_data['Repro'] , 'Resis':Tad_data['Resis'] , 'Hyper':Tad_data['Over'] , 'Others':Tad_data['Others']}
tad = {}  

H3K9mes_tads = {}

for cl in tads:
    tad[cl] = []
    H3K9mes_tads[cl] = []
    for g in chroms:
        tmp_tad = tads[cl][tads[cl]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_tad:
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() / ((end - start))
            data3 = tmp3[start:end].sum() / ((end - start))
            
            H3K9mes_tads[cl].append(data1/data3)
            if data1 > data3:
                tad[cl].append(i)
            
fig = Box_plot_4cellline([H3K9mes_tads['Repro'] , H3K9mes_tads['Resis'] , H3K9mes_tads['Hyper'] , H3K9mes_tads['Others']])
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/Clustered_TAD_Chip_intensity.pdf')   
        

p = []      
for c in ['Repro' , 'Resis' , 'Hyper' , 'Others']:
    print (c,len(tad[c]) / len(tads[c])   )
    p.append(len(tad[c]) / len(tads[c]))
    
   
    





Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (14, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.bar([1,2,3,4] , p , color = 'steelblue')
ax.text(0.6 , p[0] + 0.05 , np.round(p[0] , 7))
ax.text(2 , p[1] + 0.05 , np.round(p[1] , 7))
ax.text(3 , p[2] + 0.05 , np.round(p[2] , 7))
ax.text(3.6 , p[3] + 0.05 , np.round(p[3] , 7))


ax.set_xlim((0,5))
ax.set_ylim((0,1))
ax.set_xticks([1,2,3,4])
ax.set_xticklabels(['Repro' , 'Resis' , 'Hyper' , 'Others'] ,fontsize = 20)
ax.set_ylabel('Percentage' , fontsize = 30)
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/H3K9me3_marked_TADs_cluster_percentage_barplot.pdf')





##-----------------------Tads_classify_Chip_signal_2classifies------------------------------------   

tads = {'Repro':Tad_data['Repro'] , 'un_Repro':np.hstack((Tad_data['Resis'] , Tad_data['Over']))}
tad = {}  

H3K9mes_tads = {}

for cl in tads:
    tad[cl] = []
    H3K9mes_tads[cl] = []
    for g in chroms:
        tmp_tad = tads[cl][tads[cl]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_tad:
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() / ((end - start))
            data3 = tmp3[start:end].sum() / ((end - start))
            
            H3K9mes_tads[cl].append(data1/data3)
            if data1 > 1.4 * data3:
                tad[cl].append(i)
            
fig = Box_plot_2cellline([H3K9mes_tads['Repro'] , H3K9mes_tads['un_Repro']])
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Hitad_classify_weix/Clustered_TAD_Chip_intensity_2classifies_1.pdf')   
        














##-----------------------Tads_classify_Insulation_score------------------------------------

c1_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_CCS_noNT_fESC.txt' , dtype = tad_type)
c2_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_CCS_noNT_nofESC.txt' , dtype = tad_type)
c3_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_CCS_NT_fESC.txt' , dtype = tad_type)
c4_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_CCS_NT_nofESC.txt' , dtype = tad_type)
c5_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_noCCS_noNT_fESC.txt' , dtype = tad_type)
c6_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_noCCS_NT_fESC.txt' , dtype = tad_type)
c7_tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/Venn3_Tads/Tads_noCCS_NT_nofESC.txt' , dtype = tad_type)


tads = {'Static':c3_tad , 'Repro':np.hstack((c2_tad , c6_tad)) , 'Resis':np.hstack((c4_tad , c5_tad)) , 'Inter':np.hstack((c1_tad , c7_tad))}
     

for cl in tads:
    H3K9mes_tads = []
    for g in chroms:
        tmp_tad = tads[cl][tads[cl]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_tad:
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() / (end - start)
            data3 = tmp3[start:end].sum() / (end - start)
            if data1 > data3:
                H3K9mes_tads.append(i)
          
            
    print cl , (len(H3K9mes_tads) / len(union_tad))
    
    
    
    Insulation_score = {}
    for c in cells:
        Insulation_score[c] = []
        for i in H3K9mes_tads:
            g = i[0]
            start = i[1] // resolution
            end = i[2] // resolution
            if end - start <= 10:
                continue
            insulation1 = acquireSingleIns(HiC_Data[c][g] , start , 10 , 'divide')
            insulation2 = acquireSingleIns(HiC_Data[c][g] , end , 10 , 'divide')
            Insulation_score[c].append(insulation1)
            Insulation_score[c].append(insulation1)
    
    
    fig = Box_plot([Insulation_score['CCS'] , Insulation_score['NT5'] , Insulation_score['NT6'] , Insulation_score['F35'] , Insulation_score['F40']])
    run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/' + cl + '_marked_Tad_Insulation_Score.pdf')
    plt.close() 







CCS_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/CCS_Domain_bottom_40K_respective_stable_400K.txt' , dtype=tad_type)
CCS_Tad = CCS_Tad[CCS_Tad['chr'] != 'X']

H3K9mes_tads = [] ; n = 0
for g in chroms:
    tmp_tad = CCS_Tad[CCS_Tad['chr'] == g]
    tmp1 = Chip1[g]
    tmp3 = Input1[g]
    sum1 = tmp1.sum()
    sum3 = tmp3.sum()
    for i in tmp_tad:
        start = i[1] // 100
        end = i[2] // 100
        data1 = tmp1[start:end].sum() / (end - start)
        data3 = tmp3[start:end].sum() / (end - start)
        if data1 > data3:
            n += 1
            H3K9mes_tads.append(i)
      
        
print n , len(H3K9mes_tads) , len(CCS_Tad)
    
    
    



