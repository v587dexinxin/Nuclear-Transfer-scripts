# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 15:32:14 2018

@author: xxli
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



peaks_type = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})

motif_type = np.dtype({'names':['chr' , 'start' , 'end' , 'strand'],
                      'formats':['S8' , np.int , np.int , 'S4']})
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

def peak_index(peak):
    peak_new = {}
    for i in set(peak['chr']):
        peak_new[i] = {}
        tmp = peak[peak['chr'] == i]
        bin_size = tmp['end'].max() // 10000
        for j in range(bin_size + 1):
            peak_new[i][j] = []
        for k in tmp:
            index = k['end'] // 10000
            peak_new[i][index].append(k)
    return peak_new


def Common2_Peak(peak1 , peak2):
    peak1 = np.array(peak1 , dtype = peaks_type)
    peak2 = np.array(peak2 , dtype = peaks_type)
    peak_new = []
    for i in peak1:
        start = i['start']
        end = i['end']
        mask = (start <= peak2['end']) & (end >= peak2['start'])
        overlap = peak2[mask]
        if overlap.size != 0:
            peak_new.append(i)
        else:
            continue
    peak_new = np.array(peak_new , dtype = peaks_type)
    return peak_new


def Find_Common_Peak(peaks):
    peaks = [peak_index(peak) for peak in peaks]
    numbers = len(peaks)
    common = []
    for g in chro:
        g = 'chr' + g
        for i in peaks[0][g]:
            if peaks[0][g][i]:
                overlap = peaks[0][g][i]
                for j in range(1,numbers):
                    if i in peaks[j][g].keys():
                        overlap = Common2_Peak(overlap , peaks[j][g][i])
                    else:
                        overlap = []
                for k in overlap:
                    common.append(k)
    return common

ES1 = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\E14Tg2a.4_CTCF_NarrowPeak.bed' , dtype = peaks_type)
ES2 = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\ES_Bruce4_CTCF_NarrowPeak.bed' , dtype = peaks_type)
ES3 = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\ES_Bruce4_CTCF_NarrowPeak.bed' , dtype = peaks_type)

s_heart = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\heart_CTCF_NarrowPeak.bed' , dtype = peaks_type)
s_hindbrain = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\hindbrain_CTCF_NarrowPeak.bed' , dtype = peaks_type)
s_kidney = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\kidney_CTCF_NarrowPeak.bed' , dtype = peaks_type)
s_liver = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\liver_CTCF_NarrowPeak.bed' , dtype = peaks_type)
s_lung = np.loadtxt('G:\\data\\Chip_seq\\NarrowPeak\\lung_CTCF_NarrowPeak.bed' , dtype = peaks_type)

ctcf_motif = np.loadtxt('F:/xxli/data_analysis/motif/Ctcf_fimo/Sort_fimo.txt' ,  skiprows = 1 , usecols = (1 , 2 , 3 , 4) , dtype = motif_type)

chro = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
common = Find_Common_Peak([s_heart , s_hindbrain , s_kidney , s_liver , s_lung])

ctcf_motif = np.loadtxt('F:/xxli/data_analysis/motif/Ctcf_fimo/Sort_fimo.txt' ,  skiprows = 1 , usecols = (1 , 2 , 3 , 4) , dtype = motif_type)
chro = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
common = np.array(common , dtype = peaks_type)
motif = []
for i in ctcf_motif:
    c = i['chr']
    start = i['start']
    end = i['end']
    tmp = common[common['chr'] == c]
    mask = (start <= tmp['end']) & (end >= tmp['start'])
    overlap = tmp[mask]
    if overlap.size != 0:
        motif.append(i)
    else:
        continue



o = open('G:\\data\\Chip_seq\\NarrowPeak\\somatic_ctcf_motif.txt' , 'w')
for i in motif:
    o.writelines('\t'.join([i[0] , str(i[1]) , str(i[2]) , i[3]])  + '\n')

o.close()