# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 10:38:29 2019

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
#        if i == 1:
#            peaks[i] = Find_Common_Peak(f3,f4)
#        else:
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
            

#def Get_common2_common3_peaks(Peak1 , Peak2 , Peak3):
#    common = []
#    for i in peak1:
#        chro = i['chr']
#        start = i['start']
#        end = i['end']
#        tmp1 = peak2[peak2['chr'] == chro]
#        tmp2 = peak3[peak3['chr'] == chro]
#        mask1 = (start <= tmp1['end']) & (end >= tmp1['start'])
#        
#        overlap1 = tmp1[mask1]
#        overlap2 = tmp2[mask2]
#        if overlap1.size != 0:
#            continue
#        else:
#            
#        if (overlap1.size == 0) & (overlap2.size == 0) & (overlap3.size == 0):
#            specific.append((i['chr'] , i['start'] , i['end']))
#        else:
#            continue
#    return specific
#        
##-----------------------------------------NT_common--------------------------------------------------------------------   
datafolder = 'D:\\Workspace_New\\data\\ATAC\\peaks\\idr\\'
c1 = 'NT5' ; c2 = 'NT6'
f1 = np.loadtxt(datafolder + c1 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
f2 = np.loadtxt(datafolder + c2 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)

NT , peaks = Get_union_peaks(datafolder , ['NT5' , 'NT6'])
Write2fils_nochr('D:\\Workspace_New\\data\\ATAC\\peaks\\idr\\NT_IDR_Selected.narrowPeak' , NT)


#---------------------------------------diff_cells_common---------------------------------------------------------------

c1 = 'CCS' ; c2 = 'NT5' ; c3 = 'NT6' ; c4 = 'fESC'
CCS = np.loadtxt(datafolder + c1 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
NT5 = np.loadtxt(datafolder + c2 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
NT6 = np.loadtxt(datafolder + c3 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
fESC = np.loadtxt(datafolder + c4 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
CCS_NT5 = Find_Common2_Peak(CCS , NT5)
CCS_NT6 = Find_Common2_Peak(CCS , NT6)
CCS_fESC = Find_Common2_Peak(CCS , fESC)
NT5_fESC = Find_Common2_Peak(NT5 , fESC)
NT6_fESC = Find_Common2_Peak(NT6 , fESC)
NT5_NT6 = Find_Common2_Peak(NT5 , NT6)
CCS_NT5_fESC = Find_Common3_Peak(CCS , NT5 , fESC)
CCS_NT6_fESC = Find_Common3_Peak(CCS , NT6 , fESC)
CCS_NT5_NT6 = Find_Common3_Peak(CCS , NT5 , NT6)
NT5_NT6_fESC = Find_Common3_Peak(NT5 , NT6 , fESC)
CCS_NT5_NT6_fESC = Find_Common4_Peak(CCS , NT5 , NT6 , fESC)


union = Get_union_peaks(datafolder , ['CCS','NT5','NT6','fESC'])
CCS_specific = Get_specific(CCS , NT5 , NT6 , fESC)
NT5_specific = Get_specific(NT5 , CCS , NT6 , fESC)
NT6_specific = Get_specific(NT6 , CCS , NT5 , fESC)
fESC_specific = Get_specific(fESC , CCS , NT5 , NT6)
no_CCS = Get_specific(union[0] , CCS , CCS , CCS)
no_CCS = np.array(no_CCS , dtype = p_type)
NT5_specific = np.array(NT5_specific , dtype = p_type)
NT6_specific = np.array(NT6_specific , dtype = p_type)
fESC_specific = np.array(fESC_specific , dtype = p_type)
NT5_NT6_fESC_specific = Get_specific(no_CCS , NT5_specific , NT6_specific , fESC_specific)
  
union , peaks = Get_union_peaks(datafolder , ['CCS' , 'NT' , 'fESC'])

al = {'CCS_noNT5__noNT6_nofESC':CCS_specific , 'noCCS_NT5_noNT6_nofESC':NT5_specific , 'noCCS_noNT5_noNT6_fESC':fESC_specific ,
      'noCCS_noNT5_NT6_nofESC':NT6_specific , 'noCCS_NT5_NT6_fESC':NT5_NT6_fESC_specific}

n = 0
for i in al:
    Write2fils('D:\\Workspace_New\\data\\ATAC\\peaks\\idr\\Venn_4cells\\' + i + '.bed' , al[i])
    n += 1
    

    
print 'area1=' + str(len(CCS))
print 'area2=' + str(len(NT5))
print 'area3=' + str(len(NT6))
print 'area4=' + str(len(fESC))
print 'n12=' + str(len(CCS_NT5)) 
print 'n13=' + str(len(CCS_NT6))
print 'n14=' + str(len(CCS_fESC))
print 'n23=' + str(len(NT5_NT6))
print 'n24=' + str(len(NT5_fESC))
print 'n34=' + str(len(NT6_fESC))
print 'n123=' + str(len(CCS_NT5_NT6))
print 'n124=' + str(len(CCS_NT5_fESC))
print 'n134=' + str(len(CCS_NT6_fESC))
print 'n234=' + str(len(NT5_NT6_fESC))
print 'n1234=' + str(len(CCS_NT5_NT6_fESC))


