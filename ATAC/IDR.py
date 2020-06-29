# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 20:56:38 2018

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

cells = ['CCS' , 'NT5' , 'NT6' ,'fESC']
datatype = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['<S8' , np.int , np.int , np.float]})

itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})

def LoadPeaks(peak_file):
    '''
    '''
    itype = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})
    peaks = np.loadtxt(peak_file , dtype = itype)
    
    return peaks

def Find_Common_Peak(f1 , f2 , c):
    '''
    '''
    Common_Peak_1 = {}
    Common_Peak_2 = {}
    Peak1 = LoadPeaks(f1)
    Peak2 = LoadPeaks(f2)
    chro = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
    n = 0
    with open('D:\\Workspace_New\\data\\ATAC\\peaks\\' + c + '_R1_R2_Common_peaks.txt','w') as out:
        for g in chro:
            tmp_1 = Peak1[Peak1['chr'] == g]
            tmp_2 = Peak2[Peak2['chr'] == g]
            Common_Peak_1[g] = []
            Common_Peak_2[g] = []
            for p in tmp_1:
                mask = (tmp_2['start'] < p['end']) & (tmp_2['end'] > p['start'])
                overlap = tmp_2[mask]
                if overlap.size == 0:
                    continue
                else:
                    if overlap.size > 1:
                        n += 1
                    i = overlap[0]
                    out.writelines(p['chr'] + '\t' + str(p['start']) + '\t' + str(p['end']) + '\t' + str(overlap[0]['start'])+ \
                                   '\t' + str(overlap[0]['end']) + '\t' + str(p['score']) + '\t' + str(overlap[0]['score']) + \
                                   '\t' + str(p['p-value']) + '\t' + str(overlap[0]['p-value']) + '\n')
                    Common_Peak_1[g].append((p['chr'],p['start'],p['end'],p['name'] ,p['score'] , p['strand'] , p['value'],p['p-value'],p['q-value'],p['peak'],i['p-value']))
                    Common_Peak_2[g].append((i['chr'],i['start'],i['end'],i['name'] ,i['score'] , i['strand'] , i['value'],i['p-value'],i['q-value'],i['peak'],p['p-value']))
    with open('D:\\Workspace_New\\data\\ATAC\\peaks\\' + c + '_R1_Common_peaks.txt','w') as out:
        
        for x in Common_Peak_1.keys():
            for j in Common_Peak_1[x]:
                j = np.array(list(j),dtype = str)
                out.writelines('\t'.join(j)+'\n')
    out.close()
    with open('D:\\Workspace_New\\data\\ATAC\\peaks\\' + c + '_R2_Common_peaks.txt','w') as out:
        for x in Common_Peak_2.keys():
            for j in Common_Peak_2[x]:
                j = np.array(list(j),dtype = str)
                out.writelines('\t'.join(j)+'\n')
    out.close()
    
    return Common_Peak_1 ,Common_Peak_2,n


def Recalling_IDR_peaks(f1 , Common_fil , IDR_fil , c):
    '''
    '''
    Peak1 = LoadPeaks(f1)
    IDR_type = np.dtype({'names':['P1','P2'],
                         'formats':[np.float,np.float]})
    IDR_P = np.loadtxt(IDR_fil,dtype = IDR_type,usecols=(0,1))
    f = open(Common_fil,'r')
    with open('D:\\Workspace_New\\data\\ATAC\\peaks\\idr\\' + c + '_IDR_Selected.narrowPeak','w') as out:
        for line in f:
            l = line.strip().split()
            mask = (IDR_P['P1'] == float(l[7]))&(IDR_P['P2'] == float(l[-1]))
            tmp = IDR_P[mask]
            if tmp.size != 0:
                chro = l[0]
                start = int(l[1])
                end = int(l[2])
                mask1 = (Peak1['chr'] == chro) & (Peak1['start'] == start) & (Peak1['end'] == end)
                overlap1 = Peak1[mask1]
                if overlap1.size == 0:
                    print 1
                    continue
                else:
                    out.writelines('\t'.join([str(x) for x in overlap1[0]]) + '\n')
            else:
                continue
            


 

for c in cells:
    f1 = 'D:/Workspace_New/data/ATAC/peaks/' + c + '_R1_peaks.narrowPeak'
    f2 = 'D:/Workspace_New/data/ATAC/peaks/' + c + '_R2_peaks.narrowPeak' 
    (common_peak_1 , commom_peak_2 , n)  = Find_Common_Peak(f1 , f2 , c)
    print n


for c in ['fESC']:
    f1 = 'D:/Workspace_New/data/ATAC/peaks/' + c + '_R1_Common_peaks.txt'
    Common_fil = 'D:/Workspace_New/data/ATAC/peaks/' + c + '_R1_R2_Common_peaks.txt'
    IDR_fil = 'D:/Workspace_New/data/ATAC/peaks/' + c + '_R1_R2_selected.txt'
    Recalling_IDR_peaks(f1 , Common_fil , IDR_fil , c)

