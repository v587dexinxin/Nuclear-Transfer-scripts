# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 22:58:09 2019

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

#Nuclear transfer                  
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']

#IPSC
MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_Compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_Compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_Compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']



chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000

NT = {'CCS' : [0 , 0] , 'NT5' : [0 , 0] , 'NT6' : [0 , 0] , 'fESC' : [0 , 0]} 
IPS = {'MEF' : [0 , 0] , 'IPS_P3' : [0 , 0] , 'E14' : [0 , 0]}

for i in range(len(CCS)):
    if CCS[i]['pc'] > 0:
        NT['CCS'][0] += 1
    elif CCS[i]['pc'] < 0:
        NT['CCS'][1] += 1
    else:
        pass
    
    if NT5[i]['pc'] > 0:
        NT['NT5'][0] += 1
    elif NT5[i]['pc'] < 0:
        NT['NT5'][1] += 1
    else:
        pass
    
    if NT6[i]['pc'] > 0:
        NT['NT6'][0] += 1
    elif NT6[i]['pc'] < 0:
        NT['NT6'][1] += 1
    else:
        pass
    
    if fESC[i]['pc'] > 0:
        NT['fESC'][0] += 1
    elif fESC[i]['pc'] < 0:
        NT['fESC'][1] += 1
    else:
        pass
    
    if MEF[i]['pc'] > 0:
        IPS['MEF'][0] += 1
    elif MEF[i]['pc'] < 0:
        IPS['MEF'][1] += 1
    else:
        pass
    
    if IPS_P3[i]['pc'] > 0:
        IPS['IPS_P3'][0] += 1
    elif IPS_P3[i]['pc'] < 0:
        IPS['IPS_P3'][1] += 1
    else:
        pass
    
    if E14[i]['pc'] > 0:
        IPS['E14'][0] += 1
    elif E14[i]['pc'] < 0:
        IPS['E14'][1] += 1
    else:
        pass
    
    
x = range(4)
yA = [NT[c][0] / (NT[c][0] + NT[c][1]) for c in ['CCS' , 'NT5' , 'NT6' , 'fESC']]   
yB = [1 , 1 , 1 , 1]    
    
fig1 = Bar_plot_NT(x , yA , yB)    

run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_compartment_AB_percentage.pdf')

x = range(3)
yA = [IPS[c][0] / (IPS[c][0] + IPS[c][1]) for c in ['MEF' , 'IPS_P3' , 'E14']]   
yB = [1 , 1 , 1]    
    
fig2 = Bar_plot_IPS(x , yA , yB) 
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPS_compartment_AB_percentage.pdf')
