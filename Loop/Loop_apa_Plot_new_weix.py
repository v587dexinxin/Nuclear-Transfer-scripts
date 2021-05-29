# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 09:45:11 2021

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.special import ndtr
import math
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')


    
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
    maxi = upperpart.mean() * 5
    ## APA score
    score = avg[w,w] / lowerpart.mean()
    ## z-score
    z = (avg[w,w] - lowerpart.mean()) / lowerpart.std()
    p = 1 - ndtr(z)
    
    return avg, score, z, p, maxi
    

def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new
                
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
   
    
def Correct_VC(X, alpha):
    x = np.array(X,float)
    s1 = np.sum(x, axis = 1)
    s1 = s1 ** alpha
#    s /= np.mean(s[s!=0])
    s1[s1 == 0] = 1
    s2 = np.sum(x, axis = 0)
    s2 = s2 ** alpha
#    s2 /= np.mean(s2[s2 != 0])
    s2[s2 == 0] = 1
    return x / (s2[None, :] * s1[:, None])


def Normal_VC_Correct(NPZ):
    """
    """
    Raw_Lib = np.load(NPZ)
    Nor_Lib = {}
    for c in Raw_Lib.keys():
        Raw_M = Raw_Lib[c]
        Nor_M = Correct_VC(Raw_M, 2/3)
        
        #Recall
        Re_factor = Raw_M.mean() / Nor_M.mean()
        Nor_M = Nor_M * Re_factor
        Nor_Lib[c] = Nor_M
    
    return Nor_Lib


chros = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
res = 20000
initial_distance = 2 * math.sqrt(2)

data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})

cells = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
#HiC Data Process
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

HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'F35':F35_Lib,
            'F40':F40_Lib}







#Loop_data
CCS_speci_Loops = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/CCS_Speci_Loops_NT_fESC_merged_20K.txt' , 
                             dtype = data_type , usecols = (0 , 1 , 2))
NT_fESC_speci_Loops = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/Selected_NT_fESC_specific_loops_20K.txt' , 
                             dtype = data_type , usecols = (0 , 1 , 2))

Loops = {'CCS_speci':CCS_speci_Loops , 'NT_fESC_speci':NT_fESC_speci_Loops}


Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6                    

for cl in ['CCS_speci' , 'NT_fESC_speci']:
    for c in cells:
        P_Loops = Loops[cl]
        apa = []                    
        for g in chros:
            peaks = P_Loops[P_Loops['chr'] == g]
            pos = []
            M = HiC_Data[c][g]
            for p in peaks:
                x, y = p[1], p[2]
                if abs(y-x) < 15 * res:
                    continue
                s_l = range(p[1]//res, int(np.ceil((p[1]+20000)/float(res))))
                e_l = range(p[2]//res, int(np.ceil((p[2]+20000)/float(res))))
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
            apa.extend(tmp)
        apa = np.r_[apa]
        avg,score,z,p_value,maxi = apa_analysis(apa)
        
        fig = plt.figure(figsize = (12,12))
        plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(avg, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(avg), 0, len(avg)), vmin = 0 , vmax = 3, origin = 'lower')
        ax.set_xticks([0 , 5 , 10])
        ax.set_xticklabels(['-100K' , 'Loop'  , '100K'])
        ax.set_yticks([0 , 5 , 10])
        ax.set_yticklabels(['-100K' , 'Loop'  , '100K'])
        ax.set_title('APA score = {0:.3g}, p-value = {1:.3g} , loop_num: {2:.0f}'.format(score, p_value, len(P_Loops)))
        ax = fig.add_axes([Left + width + 0.03 , HB , 0.035 , 0.1])
        cbar = fig.colorbar(sc,cax = ax, orientation='vertical')
        cbar.set_ticks([0 , 3])
        
        run_Plot(fig , '/public/home/lixinxin/data/BDF1/HiC/Plot/Plot_new/Loop_apa/HiCCUPS/classified_by_apa_weix/F_' + c + '_' + cl + '_Loops_apa_20K_union.pdf')
        plt.close()