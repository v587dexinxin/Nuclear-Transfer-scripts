# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 21:28:32 2021

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
from scipy.special import ndtr    
import os
import math
    


def Get_loops(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    loops = []
    fileSource = os.path.join(LoopFil)
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loop_type , skiprows=1)
    Loop = Loop[Loop['chr'] != 'X']
    
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append((i['chr'].lstrip('chr') , i['start'] , i['end']))
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops



    
def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new





def Box_plot(data , cl):                
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
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'orange','linewidth':2},
            medianprops={'color':'orange','linewidth':2},
            capprops={'color':'orange','linewidth':2},
            whiskerprops={'color':'orange','linewidth':2})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'green','linewidth':2},
            medianprops={'color':'green','linewidth':2},
            capprops={'color':'green','linewidth':2},
            whiskerprops={'color':'green','linewidth':2})
    ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':2},
            medianprops={'color':'fuchsia','linewidth':2},
            capprops={'color':'fuchsia','linewidth':2},
            whiskerprops={'color':'fuchsia','linewidth':2})
    
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
    ax.set_title(cl + ',TAD_numbers:' + str(len(loops[cl])))
    ax.set_ylim((-0.25 , 8))
    
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
    ax.set_ylim((-0.25 , 7))
    
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

def apa_submatrix(M, pos, w=5):
    
    Len = M.shape[0]

    apa = []
    for i, j in pos:
        if (i-w>=0) and (i+w+1<=Len) and (j-w>=0) and (j+w+1<=Len):
            tmp = M[i-w:i+w+1, j-w:j+w+1]
            if tmp.mean()==0:
                continue
            mask = np.isnan(tmp)
            if mask.sum() > 0:
                continue
            tmp = tmp / tmp.mean()
            apa.append(tmp)
    
    return apa

def apa_analysis(apa, w=5, cw=3):
    
    avg = apa.mean(axis=0)
    lowerpart = avg[-cw:,:cw]
    upperpart = avg[:cw,-cw:]
    m = lowerpart.mean()
    maxi = upperpart.mean() * 5
    if m == 0:
        m = 0.00001
    ## APA score
    score = avg[w,w] / m
    ## z-score
    z = (avg[w,w] - lowerpart.mean()) / lowerpart.std()
    p = 1 - ndtr(z)
    
    return avg, score, z, p, maxi
    
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
    
    
    
    
    
    

loop_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                       'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})
union_type = np.dtype({'names':['chr' , 'start' , 'end' , 'CCS' , 'NTs' , 'fESC'],
                       'formats':['S8' , np.int , np.int , np.float , np.float , np.float]})


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 20000
cells = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
initial_distance = math.sqrt(2) + 0.05

##-------------------------------loop_data--------------------------

union_loop = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/union_loops_NT_fESC_merged.txt' , dtype = loop_type)

CCS_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_CCS_loops_20K_3_0.0001.txt')
NT5_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_NT5_loops_20K_3_0.0001.txt')
NT6_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_NT6_loops_20K_3_0.0001.txt')
F35_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_F35_loops_20K_3_0.0001.txt')
F40_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_F40_loops_20K_3_0.0001.txt')
NTs_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_NTs_loops_20K_3_0.000001.txt')
fESC_loop = Get_loops('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_fESC_loops_20K_3_0.000001.txt')


# NT5_1 , NT6_1 , common_NT , common_NT_2 = Get_diff_loops(NT5_loop , NT6_loop)
# F35_1 , F40_1 , common_fESC , common_fESC_2 = Get_diff_loops(F35_loop , F40_loop)

# NTs_loop = common_NT
# fESC_loop = common_fESC


#----------------------------HiC Data Process--------------------------------------
CCS_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/CCS_matrix/Cor_CCS_20K.npz' , allow_pickle=True)
CCS_Lib = Get_nan_zero_Matrix(CCS_Lib)
NT5_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT5_matrix/Cor_NT5_20K.npz' , allow_pickle=True)
NT5_Lib = Get_nan_zero_Matrix(NT5_Lib)
NT6_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NT6_matrix/Cor_NT6_20K.npz' , allow_pickle=True)
NT6_Lib = Get_nan_zero_Matrix(NT6_Lib)
F35_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F35_matrix/Cor_F35_20K.npz' , allow_pickle=True)
F35_Lib = Get_nan_zero_Matrix(F35_Lib)
F40_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/F40_matrix/Cor_F40_20K.npz' , allow_pickle=True)
F40_Lib = Get_nan_zero_Matrix(F40_Lib)

fESC_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/fESC_matrix/Cor_fESC_20K.npz' , allow_pickle=True)
fESC_Lib = Get_nan_zero_Matrix(fESC_Lib)
NTs_Lib = np.load('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/cor/NTs_matrix/Cor_NTs_20K.npz' , allow_pickle=True)
NTs_Lib = Get_nan_zero_Matrix(NTs_Lib)


HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'F35':F35_Lib,
            'F40':F40_Lib,
            'NTs':NTs_Lib,
            'fESC':fESC_Lib}




##--------------------------Loops_classify_new--------------------------------




union_apa = []
for loop in union_loop:
    g = loop['chr']
    x, y = loop[1], loop[2]
    apa_scores = []
    for c in ['CCS' , 'NTs' , 'fESC']:
        M = HiC_Data[c][g]
        matrix = M[x // res : y // res , x // res : y // res]
        if len(np.nonzero(matrix)[0]) / (len(matrix)**2)  <= 0.3:
            print c , loop
            continue
        pos = []
        apa = []
        s_l = range(x//res, int(np.ceil((x+20000)/float(res))))
        e_l = range(y//res, int(np.ceil((y+20000)/float(res))))
        si, ei = None, None
        for st in s_l:
            for et in e_l:
                if (st < M.shape[0]) and (et < M.shape[0]):
                    if si is None:
                        si, ei = st, et
                    else:
                        if M[st,et]>M[si,ei]:
                            si,ei = st,et
        if not si is None:
            if si < ei:
                pos.append((si,ei))
            else:
                pos.append((ei,si))
        tmp = apa_submatrix(M, pos)
        if len(tmp) == 0:
            score = 0
        else:
            
            apa.extend(tmp)
            apa = np.r_[apa]
            avg,score,z,p,maxi = apa_analysis(apa)
        apa_scores.append(score)
    if len(apa_scores) == 3:
        if not (np.isnan(apa_scores[0]) or np.isnan(apa_scores[1]) or np.isnan(apa_scores[2]) or np.isinf(apa_scores[0]) or np.isinf(apa_scores[1]) or np.isinf(apa_scores[2])):
        # if (((apa_scores[1] - apa_scores[2]) > -1) and ((apa_scores[1] - apa_scores[2]) < 1)) and (((apa_scores[3] - apa_scores[4]) > -1) and ((apa_scores[3] - apa_scores[4]) < 1)):
            union_apa.append((g , x , y , apa_scores[0] , apa_scores[1] , apa_scores[2]))
union_apa = np.array(union_apa , dtype = union_type)
        
 
        
 
##Specific
CCS_speci = []
for i in union_apa:
    g = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = CCS_loop[CCS_loop['chr'] == g]
    tmp2 = fESC_loop[fESC_loop['chr'] == g]
    tmp3 = NTs_loop[NTs_loop['chr'] == g]
    
    mask1 = (tmp1['start'] >= (start - res)) & (tmp1['start'] <= (start + res)) & (tmp1['end'] >= (end - res)) & (tmp1['end'] <= (end + res))
    mask2 = (tmp2['start'] >= (start - res)) & (tmp2['start'] <= (start + res)) & (tmp2['end'] >= (end - res)) & (tmp2['end'] <= (end + res))
    mask3 = (tmp3['start'] >= (start - res)) & (tmp3['start'] <= (start + res)) & (tmp3['end'] >= (end - res)) & (tmp3['end'] <= (end + res))
    
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]
    overlap3 = tmp3[mask3]
    
    if (overlap1.size != 0) and (overlap2.size == 0) and (overlap3.size == 0):
        CCS_speci.append(i)
CCS_speci = np.array(CCS_speci , dtype = union_apa.dtype)

fESC_speci = []
for i in union_apa:
    g = i['chr']
    start = i['start']
    end = i['end']
    
    tmp1 = CCS_loop[CCS_loop['chr'] == g]
    tmp2 = fESC_loop[fESC_loop['chr'] == g]
    tmp3 = NTs_loop[NTs_loop['chr'] == g]
    
    mask1 = (tmp1['start'] >= (start - res)) & (tmp1['start'] <= (start + res)) & (tmp1['end'] >= (end - res)) & (tmp1['end'] <= (end + res))
    mask2 = (tmp2['start'] >= (start - res)) & (tmp2['start'] <= (start + res)) & (tmp2['end'] >= (end - res)) & (tmp2['end'] <= (end + res))
    mask3 = (tmp3['start'] >= (start - res)) & (tmp3['start'] <= (start + res)) & (tmp3['end'] >= (end - res)) & (tmp3['end'] <= (end + res))
    
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]
    overlap3 = tmp3[mask3]
    
    if (overlap1.size == 0) and (overlap2.size != 0) and (overlap3.size != 0):
        fESC_speci.append(i)
fESC_speci = np.array(fESC_speci , dtype = union_apa.dtype)





## Classify
A_B_Repro = [] ; A_B_Resis = [] ; A_B_Over = []
B_A_Repro = [] ; B_A_Resis = [] ; B_A_Over = [] ; Others = []


for i in union_apa:
    g = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = CCS_loop[CCS_loop['chr'] == g]
    tmp2 = fESC_loop[fESC_loop['chr'] == g]
    mask1 = (tmp1['start'] >= (start - res)) & (tmp1['start'] <= (start + res)) & (tmp1['end'] >= (end - res)) & (tmp1['end'] <= (end + res))
    mask2 = (tmp2['start'] >= (start - res)) & (tmp2['start'] <= (start + res)) & (tmp2['end'] >= (end - res)) & (tmp2['end'] <= (end + res))
    overlap1= tmp1[mask1]
    overlap2 = tmp2[mask2]

    if ((i['CCS'] - i['fESC']) >= -0.3) and ((i['CCS'] - i['fESC']) <= 0.3):
        Others.append(i)
    elif ((overlap1.size != 0) and (overlap2.size == 0) and (i['CCS'] - i['fESC'] > 0.3 )):
        if  ((i['NTs'] - i['CCS']) == 0 ) or ((i['fESC'] - i['CCS']) == 0):
            print i
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 1.2:
            A_B_Over.append(i)
        elif ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 0.7) and ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 1.2):
            A_B_Repro.append(i)
        elif (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 0.7:
            A_B_Resis.append(i)
    elif (overlap1.size == 0) and (overlap2.size != 0) and (i['CCS'] - i['fESC'] < -0.3 ):
        if  ((i['NTs'] - i['CCS']) == 0 ) or ((i['fESC'] - i['CCS']) == 0):
            print i
        if (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 1.2:
            B_A_Over.append(i)
        elif ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) > 0.7) and ((i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 1.2):
            B_A_Repro.append(i)
        elif (i['NTs'] - i['CCS']) / (i['fESC'] - i['CCS']) <= 0.7:
            B_A_Resis.append(i)





A_B_Repro = np.array(A_B_Repro , dtype = union_apa.dtype)
A_B_Resis = np.array(A_B_Resis , dtype = union_apa.dtype)
A_B_Over = np.array(A_B_Over , dtype = union_apa.dtype)

B_A_Repro = np.array(B_A_Repro , dtype = union_apa.dtype)
B_A_Resis = np.array(B_A_Resis , dtype = union_apa.dtype)
B_A_Over = np.array(B_A_Over , dtype = union_apa.dtype)
Others = np.array(Others , dtype = union_apa.dtype)



loops = {'A_B_Repro':A_B_Repro , 'A_B_Resis':A_B_Resis , 'A_B_Over':A_B_Over ,
        'B_A_Repro':B_A_Repro , 'B_A_Resis':B_A_Resis , 'B_A_Over':B_A_Over ,
        'Others':Others, 
        'Repro':np.hstack((A_B_Repro , B_A_Repro)) , 'Resis':np.hstack((A_B_Resis , B_A_Resis)) , 
        'Over':np.hstack((A_B_Over , B_A_Over)) , 
        'CCS_Speci':CCS_speci ,'fESC_Speci':fESC_speci}


##---------------------------plot_3celllines----------------------------------------     

for cl in loops:
    
    out = open('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/' + cl + '_Loops_NT_fESC_merged_20K.txt' , 'w')
    for i in loops[cl]:
        out.writelines('\t'.join([i['chr'] , str(i['start']) , str(i['end']) ]) + '\n')
    out.close()
    
    
    fig = Box_plot_3cellline([loops[cl]['CCS'] , loops[cl]['NTs'] , loops[cl]['fESC']] , cl)
    run_Plot(fig , '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/plot/' + cl + '_Loops_apa_boxplot_3cellines.pdf')
    
    print cl , len(loops[cl])

    plt.close()
    

  
##---------------------------plot_3celllines_merged---------------------------------------
fig = Box_plot_3cellline_merged([loops['Repro']['CCS'] , loops['Repro']['NTs'] , loops['Repro']['fESC'] , loops['Resis']['CCS'] , loops['Resis']['NTs'] , loops['Resis']['fESC'] , loops['Over']['CCS'] , loops['Over']['NTs'] , loops['Over']['fESC'] , loops['Others']['CCS'] , loops['Others']['NTs'] , loops['Others']['fESC']])
run_Plot(fig , '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs_new/classify_loops/plot/Loop_apa_3celllines_allclassify.pdf')
plt.close()


