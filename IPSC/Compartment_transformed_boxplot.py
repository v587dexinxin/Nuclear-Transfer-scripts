# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 15:33:14 2019

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
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['skyblue' , 'k' , 'yellow'])
#my_cmap.set_bad('#2672a1')

pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['MEF' , 'IPS_P3' , 'E14'] , 
                    'formats':[np.float , np.float , np.float , np.float]})
 
def PC_changed(MEF , IPS_P3 , E14):
    A_B = []
    B_A = []
    for i in range(len(MEF)):
        if MEF[i]['chr'] in chroms:
            if (MEF[i]['pc'] > 0) and (E14[i]['pc'] < 0):
                A_B.append((MEF[i]['pc'] , IPS_P3[i]['pc'] , E14[i]['pc']))
            elif (MEF[i]['pc'] < 0) and (E14[i]['pc'] > 0):
                B_A.append((MEF[i]['pc'] , IPS_P3[i]['pc'] , E14[i]['pc']))
    A_B = np.array(A_B , dtype = data_type)
    B_A = np.array(B_A , dtype = data_type)
    return A_B , B_A

def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'green','linewidth':2},
            medianprops={'color':'green','linewidth':2},
            capprops={'color':'green','linewidth':2},
            whiskerprops={'color':'green','linewidth':2})

    ax.plot([0.5,3.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels(['MEF' , 'IPS_P3' , 'E14'] , fontsize = 28)
    ax.set_xlim((0.5 , 3.5))
    ax.set_ylim((-0.1 , 0.1))
    return fig
             
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
                   
                    
MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_Compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

A_B , B_A = PC_changed(MEF , IPS_P3 , E14)

data1 = [A_B['MEF'] , A_B['IPS_P3'] , A_B['E14']]
data2 = [B_A['MEF'] , B_A['IPS_P3'] , B_A['E14']]

fig1 = Box_plot(data1)
fig2 = Box_plot(data2)
run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_compartment_A_B_200K.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_compartment_B_A_200K.pdf')