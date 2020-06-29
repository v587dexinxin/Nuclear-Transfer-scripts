# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 09:53:24 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles


tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()
    
def Load_Tads(TadFil):
    Tads = np.loadtxt(TadFil , dtype = tad_type)
    return Tads

def union_Tads(data_list):
    Tads_all = []
    for c in data_list:
        for i in c:
            Tads_all.append(i)
    Tads_all.sort(key = lambda x:(x[0],x[1],x[2]))
    Tads_all = np.array(Tads_all ,dtype = tad_type)
    union_tads = []
    for g in chrom:
        Tads = Tads_all[Tads_all['chr'] == g]
        Tads_number = len(Tads)
        for i in range(Tads_number-1 , -1 , -1):
            start = Tads[i]['start']
            prev_ss = Tads[i-1]['start'] - 40000
            prev_se = Tads[i-1]['start'] + 40000
            if (start >= prev_ss) and (start <= prev_se):
                Tads = list(Tads)
                Tads.remove(Tads[i])
                Tads = np.array(Tads,dtype = tad_type)
            else:
                pass
        union_tads.extend(list(Tads))
            
    union_tads = np.array(union_tads , dtype = tad_type) 
    return union_tads
    
    
    
def Common_Tads(data_list , union_tads):
    commom = []
    for i in union_tads:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        n = 0
        for j in data_list:
            tmp = j[j['chr'] == chro]
            mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
            overlap = tmp[mask]
            if overlap.size != 0:
                n += 1
            else:
                pass
        if n == len(data_list):
            commom.append(i)
    commom = np.array(commom , dtype = tad_type)
    return commom

def Common2_Tads(tads1 , tads2):
    commom = []
    n1 = 0
    n2 = 0
    for i in tads1:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tads2[tads2['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size != 0:
            commom.append(i)
            n1 += 1
    for i in tads2:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tads1[tads1['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size != 0:
            n2 += 1
    commom = np.array(commom , dtype = tad_type)
    return commom , n1 , n2
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
            
            
TadFolder = 'H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19']
cells = ['NT5','NT6','fESC']

CCS = Load_Tads(os.path.join(TadFolder , 'CCS_40K_allreps_union_small_domain.txt'))
NT5 = Load_Tads(os.path.join(TadFolder , 'NT5_40K_allreps_union_small_domain.txt'))
NT6 = Load_Tads(os.path.join(TadFolder , 'NT6_40K_allreps_union_small_domain.txt'))
fESC = Load_Tads(os.path.join(TadFolder , 'fESC_40K_allreps_union_small_domain.txt'))

data_list = [NT5 , NT6 , fESC]
NT_ESC_union = union_Tads(data_list)

CCS_NTESCunion_common , n1 , n2 = Common2_Tads(CCS , NT_ESC_union)

NT_ESC_common = Common_Tads(data_list , union_Tads(data_list))

CCS_NTESC_common , n3 , n4 = Common2_Tads(CCS , NT_ESC_common)


fig = plt.figure(figsize = (10, 10))
venn2(subsets=(len(CCS) - n3 , len(NT_ESC_common) - n4 , len(CCS_NTESC_common)), set_labels=('CCS', 'NT_fESC_common'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\CCS_NTfESC(common)_Venn2.pdf')


fig = plt.figure(figsize = (10, 10))
venn2(subsets=(len(CCS) - n1 , len(NT_ESC_union) - n2 , len(CCS_NTESCunion_common)), set_labels=('CCS', 'NT_fESC_union'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\CCS_NTfESC(union)_Venn2.pdf')



  