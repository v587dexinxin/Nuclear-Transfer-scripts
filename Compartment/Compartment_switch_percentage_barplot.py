# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 16:13:49 2020

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
data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})

def Bar_plot_NT(x ,y):
    left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (8, 12))
    ax = fig.add_axes(size_axes)
    ax.bar(x , y , color = 'skyblue')
    xticks = x
    labels = ['A-B' , 'B-A']
    for a, b in zip(x, y):
        plt.text(a, b + 0.01, '%.4f' % b, ha='center', va='bottom', fontsize=20)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('compartment type' , fontsize = 20 )
    ax.set_xlim((-0.5, 1.5))
    ax.set_ylim((0 , 0.2))
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


c1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c7 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c8 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)


x = range(2)

a = (len(c1) + len(c2) + len(c3) + len(c4)) / len(CCS)
b = (len(c5) + len(c6) + len(c7) + len(c8)) / len(CCS)

y = [a , b]

fig = Bar_plot_NT(x , y)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\Percentage_compartment_switch.pdf')



