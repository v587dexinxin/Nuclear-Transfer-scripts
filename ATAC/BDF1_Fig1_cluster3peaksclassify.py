# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 09:53:02 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import csv
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice
import os

cells = ['CCS' , 'NT5' , 'NT6' ,'fESC']

p_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['<S8' , np.int , np.int]})

datatype = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['<S8' , np.int , np.int , np.float]})

peaks_type = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})

def Find_Common2_Peak(Peak1,Peak2):
    '''
    '''
    Common_Peak = []

    for i in Peak1:
        chro = i['chr']
        tmp_2 = Peak2[Peak2['chr'] == chro]
        mask2 = (i['start'] <= tmp_2['end']) & (i['end'] >= tmp_2['start'])
        overlap2 = tmp_2[mask2]
        if overlap2.size == 0:
            continue
        else:
            Common_Peak.append(i)
                
    Common_Peak = np.array(Common_Peak , dtype = peaks_type)    
    return Common_Peak


def Find_Common3_Peak(Peak1,Peak2,Peak3):
    '''
    '''
    Common_Peak = []

    for i in Peak1:
        chro = i['chr']
        tmp_2 = Peak2[Peak2['chr'] == chro]
        tmp_3 = Peak3[Peak3['chr'] == chro]
        mask2 = (i['start'] <= tmp_2['end']) & (i['end'] >= tmp_2['start'])
        overlap2 = tmp_2[mask2]
        if overlap2.size == 0:
            continue
        else:
            mask3 = (i['start'] <= tmp_3['end']) & (i['end'] >= tmp_3['start'])
            overlap3 = tmp_3[mask3]
            mask4 = (overlap2[0]['start'] <= tmp_3['end']) & (overlap2[0]['end'] >= tmp_3['start'])
            overlap4 = tmp_3[mask4]
            if overlap3.size == 0 and overlap4.size == 0:
                continue
            else:
                Common_Peak.append(i)
                
    Common_Peak = np.array(Common_Peak , dtype = peaks_type)   
    return Common_Peak





def Find_Common4_Peak(Peak1,Peak2,Peak3,Peak4):
    '''
    '''
    Common_Peak = []

    for i in Peak1:
        chro = i['chr']
        tmp_2 = Peak2[Peak2['chr'] == chro]
        tmp_3 = Peak3[Peak3['chr'] == chro]
        tmp_4 = Peak4[Peak4['chr'] == chro]
        mask2 = (i['start'] <= tmp_2['end']) & (i['end'] >= tmp_2['start'])
        overlap2 = tmp_2[mask2]
        if overlap2.size == 0:
            continue
        else:
            mask3 = (i['start'] <= tmp_3['end']) & (i['end'] >= tmp_3['start'])
            overlap3 = tmp_3[mask3]
            if overlap3.size == 0:
                continue
            else:
                mask4 = (i['start'] <= tmp_4['end']) & (i['end'] >= tmp_4['start'])
                overlap4 = tmp_4[mask4]
                if overlap4.size == 0:
                    continue
                else:
                    Common_Peak.append(i)
                
    Common_Peak = np.array(Common_Peak , dtype = peaks_type)   
    return Common_Peak


def Write2fils(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('chr' + i[0] + '\t' + '\t'.join(i[1:])+'\n')
    out.close()
    
    
def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()
    
    
    
def Get_specific(peak1 , peak2 , peak3 , peak4):
    specific = []
    for i in peak1:
        chro = i['chr']
        start = i['start']
        end = i['end']
        tmp1 = peak2[peak2['chr'] == chro]
        tmp2 = peak3[peak3['chr'] == chro]
        tmp3 = peak4[peak4['chr'] == chro]
        mask1 = (start <= tmp1['end']) & (end >= tmp1['start'])
        mask2 = (start <= tmp2['end']) & (end >= tmp2['start'])
        mask3 = (start <= tmp3['end']) & (end >= tmp3['start'])
        overlap1 = tmp1[mask1]
        overlap2 = tmp2[mask2]
        overlap3 = tmp3[mask3]
        if (overlap1.size == 0) & (overlap2.size == 0) & (overlap3.size == 0):
            specific.append((i['chr'] , i['start'] , i['end']))
        else:
            continue
    return specific

def Get_union_peaks(datafolder , cell):
    """
    cell : List ,The List of the cells'name , e.g.['fESC' , 'NT5NT6' , 'CCS']
    Len : The window of loop_overlap
    """
    peaks = {}
    union = np.loadtxt(datafolder + cell[0] + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
    peaks[0] = np.loadtxt(datafolder + cell[0] + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
    for i in range(1,len(cell)):
        peaks[i] = np.loadtxt(datafolder + cell[i] + '_IDR_Selected.narrowPeak' , dtype = peaks_type)           
        for j in peaks[i]:
            mask = (j['start'] <=  union[union['chr'] == j['chr']]['end']) & (j['end'] >=  union[union['chr'] == j['chr']]['start'])
            tmp = union[union['chr'] == j['chr']][mask]
            if tmp.size == 0:
                t = list(union)
                t.append(j)
                union = np.array(t , dtype = peaks_type)
            else:
                continue
    return union,peaks
    
##-----------------------------------------NT_union--------------------------------------------------------------------   
datafolder = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\'
#c1 = 'NT5' ; c2 = 'NT6'
#f1 = np.loadtxt(datafolder + c1 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
#f2 = np.loadtxt(datafolder + c2 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
NT , peaks = Get_union_peaks(datafolder , ['NT5' , 'NT6'])
NT_peaks = []
for i in NT:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp1 = peaks[0][peaks[0]['chr'] == chro]
    tmp2 = peaks[1][peaks[1]['chr'] == chro]
    mask1 = (tmp1['start'] <= end) & (tmp1['end'] >= start)
    mask2 = (tmp2['start'] <= end) & (tmp2['end'] >= start)
    overlap1 = tmp1[mask1]
    overlap2 = tmp2[mask2]
    if (overlap1.size != 0) and (overlap2.size != 0):
        score = (overlap1['score'].mean() + overlap2['score'].mean()) / 2
    elif (overlap1.size != 0) and (overlap2.size == 0):
        score = overlap1['score'].mean()
    elif (overlap1.size == 0) and (overlap2.size != 0):
        score = overlap2['score'].mean()
    else:
        print i
    NT_peaks.append((chro , start , end , score))
NT_peaks = np.array(NT_peaks , dtype = datatype)
##fESC_peaks
fESC_peaks = np.loadtxt(os.path.join(datafolder , 'fESC_IDR_Selected.narrowPeak'), dtype = datatype , usecols = (0 , 1 , 2 , 4))


##--------------------------------------cluster3_peaks----------------------------------------------------
cluster3 = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster3.bed' , dtype = p_type)
c1 = [] ; c2 = [] ; c3 = []
for i in cluster3:
    chro = i['chr'].lstrip('chr')
    start = i['start']
    end = i['end']
    tmp1 = NT_peaks[NT_peaks['chr'] == chro]
    tmp2 = fESC_peaks[fESC_peaks['chr'] == chro]
    mask1 = (tmp1['start'] <= end) & (tmp1['end'] >= start)
    mask2 = (tmp2['start'] <= end) & (tmp2['end'] >= start)
    overlap1 = tmp1[mask1]
    overlap2 = tmp2[mask2]
    if (overlap1.size != 0) and (overlap2.size != 0):
        s_NT = overlap1['score'].mean()
        s_fESC = overlap2['score'].mean()
        if s_fESC > s_NT * 10:
            c3.append(i)
        elif s_NT > s_fESC * 10:
            c2.append(i)
        else:
            c1.append(i)
    elif (overlap1.size != 0) and (overlap2.size == 0):
        c2.append(i)
    elif (overlap1.size == 0) and (overlap2.size != 0):
        c3.append(i)
    else:
        c1.append(i)
            
w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster3_c1-1.5.bed' , 'w')
for i in c1:
    w.writelines('\t'.join([str(x) for x in i]) + '\n')
w.close()
            