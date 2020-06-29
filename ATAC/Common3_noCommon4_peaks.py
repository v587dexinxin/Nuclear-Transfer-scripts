# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 16:49:49 2019

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





def Find_Common3_noCommon4_Peak(Peak1,Peak2,Peak3,Peak4):
    '''
    '''
    common = Find_Common4_Peak(Peak1)
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
            mask4 = (overlap2[0]['start'] <= tmp_3['end']) & (overlap2[0]['end'] >= tmp_3['start'])
            overlap4 = tmp_3[mask4]
            if overlap3.size == 0 and overlap4.size == 0:
                continue
            else:
                Common_Peak.append(i)
                
    Common_Peak = np.array(Common_Peak , dtype = peaks_type)   
    return Common_Peak

cells = ['CCS' , 'NT5' , 'NT6' ,'fESC']

p_type = ({'names':['chr' , 'start' , 'end'],
             'formats':['<S8' , np.int , np.int]})

datatype = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['<S8' , np.int , np.int , np.float]})

peaks_type = np.dtype({'names':['chr','start','end','name' ,'score' , 'strand' , 'value','p-value','q-value','peak'],
                      'formats':['<S8' , np.int , np.int , '<S64', np.float , '<S5' , np.float, np.float, np.float , np.int]})
    
    c1 = 'CCS' ; c2 = 'NT5' ; c3 = 'NT6' ; c4 = 'fESC'
CCS = np.loadtxt(datafolder + c1 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
NT5 = np.loadtxt(datafolder + c2 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
NT6 = np.loadtxt(datafolder + c3 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)
fESC = np.loadtxt(datafolder + c4 + '_IDR_Selected.narrowPeak' , dtype = peaks_type)