# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 22:04:33 2019

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


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float]})



def get_AB_percent(data):
    data_1 = {'CCS' : [0 , 0] , 'NT5' : [0 , 0] , 'NT6' : [0 , 0] , 'fESC' : [0 , 0]}
    for i in data:
        if i[0] > 0:
            data_1['CCS'][0] += 1
        elif i[0] < 0:
            data_1['CCS'][1] += 1
        
        if i[1] > 0:
            data_1['NT5'][0] += 1
        elif i[1] < 0:
            data_1['NT5'][1] += 1
            
        if i[2] > 0:
            data_1['NT6'][0] += 1
        elif i[2] < 0:
            data_1['NT6'][1] += 1
            
        if i[3] > 0:
            data_1['fESC'][0] += 1
        elif i[3] < 0:
            data_1['fESC'][1] += 1
    A_percent = {}
    for c in data_1:
        A_percent[c] = data_1[c][0] / (data_1[c][0] + data_1[c][1])
    return A_percent
    
    
    
def get_IPSC_AB_percent(data):
    data_1 = {'MEF' : [0 , 0] , 'IPS_P3' : [0 , 0] , 'E14' : [0 , 0]}
    for i in data:
        if i[0] > 0:
            data_1['MEF'][0] += 1
        elif i[0] < 0:
            data_1['MEF'][1] += 1
        
        if i[1] > 0:
            data_1['IPS_P3'][0] += 1
        elif i[1] < 0:
            data_1['IPS_P3'][1] += 1
            
        if i[2] > 0:
            data_1['E14'][0] += 1
        elif i[2] < 0:
            data_1['E14'][1] += 1
            

    A_percent = {}
    for c in data_1:
        A_percent[c] = data_1[c][0] / (data_1[c][0] + data_1[c][1])
    return A_percent
    
def Bar_plot_NT(x , yA , yB):
    left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.bar(x , yB , color = 'darkblue')
    ax.bar(x , yA , color = 'darkorange')
    xticks = x
    labels = ['CCS' , 'NT5' , 'NT6' , 'fESC']
    for a, b in zip(x, yA):
        plt.text(a, b + 0.05, '%.4f' % b, ha='center', va='bottom', fontsize=20)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('Cell type' , fontsize = 20 )
    ax.set_xlim((-0.5, 3.5))
    ax.set_ylim((0 , 1.01))
    return fig
    
def Bar_plot_IPS(x , yA , yB):
    left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.bar(x , yB , color = 'darkblue')
    ax.bar(x , yA , color = 'darkorange')
    xticks = x
    labels = ['MEF' , 'IPS_P3' , 'E14']
    for a, b in zip(x, yA):
        plt.text(a, b + 0.05, '%.4f' % b, ha='center', va='bottom', fontsize=20)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('Cell type' , fontsize = 20 )
    ax.set_xlim((-0.5, 2.5))
    ax.set_ylim((0 , 1.01))
    return fig        
   
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

                    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']


MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


A_B = []
B_A = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_fESC[i]['pc'] < 0):
            A_B.append([tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc']])
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_fESC[i]['pc'] > 0):
            B_A.append([tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc']])

A_B_1 = get_AB_percent(A_B)
B_A_1 = get_AB_percent(B_A)


x = range(4)
yA = [B_A_1[c] for c in ['CCS' , 'NT5' , 'NT6' , 'fESC']]   
yB = [1 , 1 , 1 , 1]   

fig = Bar_plot_NT(x , yA , yB)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_compartment_B_A_percentage.pdf')


##IPSC
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            IPSC_A_B.append([tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']])
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            IPSC_B_A.append([tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']])

IPSC_A_B_1 = get_IPSC_AB_percent(IPSC_A_B)
IPSC_B_A_1 = get_IPSC_AB_percent(IPSC_B_A)

x = range(3)
yA = [IPSC_B_A_1[c] for c in ['MEF' , 'IPS_P3' , 'E14']]   
yB = [1 , 1 , 1]   

fig = Bar_plot_IPS(x , yA , yB)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_compartment_B_A_percentage.pdf')

