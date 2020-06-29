# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 19:14:40 2019

@author: han-luo
"""

from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['mediumblue' , 'white' , 'red'])
#my_cmap.set_bad('#2672a1')
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')

def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

def plot_heatmap(matrix , vmin , vmax):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    x = ['MEF','IPS_P3','E14']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 10)
    ax.set_ylabel('c1:' + str(len(IPSC_c1)) + ' , c2:' + str(len(IPSC_c2)) + ' , c3:' + str(len(IPSC_c3)) + ' , c4:' + str(len(IPSC_c4)))
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Compartment numbers:' + str(len(matrix)),fontsize = 30)
    ##Colorbar
    ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
    cbar = fig.colorbar(im,cax = ax, orientation='vertical')
    cbar.set_ticks([vmin , vmax])
    
    return fig

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['chr' , 'start' , 'end' , 'CCS_pc1' , 'NT5_pc1' , 'NT6_pc1' , 'fESC_pc1']) + '\n')
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
   
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

pc_type = np.dtype({'names':['chr' , 'score'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['chr' , 'start' , 'CCS' , 'NT5' , 'NT6' , 'fESC'] , 
                      'formats':['S8' , np.int , np.float , np.float , np.float , np.float]})
data_type_1 = np.dtype({'names':['chr' , 'start' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'] , 
                        'formats':['S8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

c1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 3 , 4 , 5 , 6) , dtype = data_type)
c2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 3 , 4 , 5 , 6) , dtype = data_type)
c3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 3 , 4 , 5 , 6) , dtype = data_type)
c4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 3 , 4 , 5 , 6) , dtype = data_type)


MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']

IPSC_c1 = [] ; IPSC_c2 = [] ; IPSC_c3 = [] ; IPSC_c4 = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    tmp_c1 = c1[c1['chr'] == g]
    tmp_c2 = c2[c2['chr'] == g]
    tmp_c3 = c3[c3['chr'] == g]
    tmp_c4 = c4[c4['chr'] == g]
    for i in tmp_c1:
        start = i['start'] // 200000
        MEF_score = tmp_MEF[start]['score']
        IPS_P3_score = tmp_IPS_P3[start]['score']
        E14_score = tmp_E14[start]['score']
        IPSC_c1.append((g , i['start'] , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , MEF_score , IPS_P3_score , E14_score))
        
    for i in tmp_c2:
        start = i['start'] // 200000
        MEF_score = tmp_MEF[start]['score']
        IPS_P3_score = tmp_IPS_P3[start]['score']
        E14_score = tmp_E14[start]['score']
        IPSC_c2.append((g , i['start'] , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , MEF_score , IPS_P3_score , E14_score))
    
    for i in tmp_c3:
        start = i['start'] // 200000
        MEF_score = tmp_MEF[start]['score']
        IPS_P3_score = tmp_IPS_P3[start]['score']
        E14_score = tmp_E14[start]['score']
        IPSC_c3.append((g , i['start'] , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , MEF_score , IPS_P3_score , E14_score))
        
    for i in tmp_c4:
        start = i['start'] // 200000
        MEF_score = tmp_MEF[start]['score']
        IPS_P3_score = tmp_IPS_P3[start]['score']
        E14_score = tmp_E14[start]['score']
        IPSC_c4.append((g , i['start'] , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , MEF_score , IPS_P3_score , E14_score))
        
IPSC_c1 = np.array(IPSC_c1 , dtype = data_type_1)
IPSC_c2 = np.array(IPSC_c2 , dtype = data_type_1)
IPSC_c3 = np.array(IPSC_c3 , dtype = data_type_1)
IPSC_c4 = np.array(IPSC_c4 , dtype = data_type_1)
    
matrix = np.zeros((len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3) + len(IPSC_c4) , 3))
for i in range(len(IPSC_c1)):
    matrix[i][0] = IPSC_c1[i]['MEF']
    matrix[i][1] = IPSC_c1[i]['IPS_P3']
    matrix[i][2] = IPSC_c1[i]['E14']

for i in range(len(IPSC_c1) , len(IPSC_c1) + len(IPSC_c2)):
    matrix[i][0] = IPSC_c2[i - len(IPSC_c1)]['MEF']
    matrix[i][1] = IPSC_c2[i - len(IPSC_c1)]['IPS_P3']
    matrix[i][2] = IPSC_c2[i - len(IPSC_c1)]['E14']
    
for i in range(len(IPSC_c1) + len(IPSC_c2) , len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3)):
    matrix[i][0] = IPSC_c3[i - (len(IPSC_c1) + len(IPSC_c2))]['MEF']
    matrix[i][1] = IPSC_c3[i - (len(IPSC_c1) + len(IPSC_c2))]['IPS_P3']
    matrix[i][2] = IPSC_c3[i - (len(IPSC_c1) + len(IPSC_c2))]['E14']


for i in range(len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3) , len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3) + len(IPSC_c4)):
    matrix[i][0] = IPSC_c4[i - (len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3))]['MEF']
    matrix[i][1] = IPSC_c4[i - (len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3))]['IPS_P3']
    matrix[i][2] = IPSC_c4[i - (len(IPSC_c1) + len(IPSC_c2) + len(IPSC_c3))]['E14']        
        
        
matrix = Z_score(matrix)
fig = plot_heatmap(matrix , -1.5 , 1.5)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\Compartment_classify_heatmap.pdf')











