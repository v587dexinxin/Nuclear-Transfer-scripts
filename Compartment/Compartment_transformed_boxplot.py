# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 20:23:47 2019

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
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float]})
 
def PC_changed(CCS , NT5 , NT6 ,fESC):
    A_B = []
    B_A = []
    for i in range(len(CCS)):
        if CCS[i]['chr'] in chroms:
            if (CCS[i]['pc'] > 0) and (fESC[i]['pc'] < 0):
                A_B.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , fESC[i]['pc']))
            elif (CCS[i]['pc'] < 0) and (fESC[i]['pc'] > 0):
                B_A.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , fESC[i]['pc']))
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
            boxprops={'color': 'orange','linewidth':2},
            medianprops={'color':'orange','linewidth':2},
            capprops={'color':'orange','linewidth':2},
            whiskerprops={'color':'orange','linewidth':2})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'green','linewidth':2},
            medianprops={'color':'green','linewidth':2},
            capprops={'color':'green','linewidth':2},
            whiskerprops={'color':'green','linewidth':2})
    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F40' ] , fontsize = 28)
    ax.set_xlim((0.5 , 4.5))
    ax.set_ylim((-0.1 , 0.1))
    return fig
             
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
                   
                    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

A_B , B_A = PC_changed(CCS , NT5 , NT6 , fESC)

data1 = [A_B['CCs'] , A_B['NT5'] , A_B['NT6'] , A_B['F40']]
data2 = [B_A['CCs'] , B_A['NT5'] , B_A['NT6'] , B_A['F40']]

fig1 = Box_plot(data1)
fig2 = Box_plot(data2)
run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\compartment_A_B_200K.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\compartment_B_A_200K.pdf')
