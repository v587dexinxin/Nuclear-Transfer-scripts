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
                    'formats':['U8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F35' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float , np.float]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']


def Load_PC(PC_fil):
    '''
    '''
    data = np.loadtxt(PC_fil , dtype = pc_type)
    data_new = []
    for g in chroms:
        tmp = data[data['chr'] == g]
        for i in tmp:
            data_new.append(i)
    data_new = np.array(data_new , dtype = data.dtype)
    
    return data_new
    
    
    
def PC_changed(CCS , NT5 , NT6 , F35 , F40):
    A_B = []
    B_A = []
    for i in range(len(CCS)):
        if CCS[i]['chr'] in chroms:
            if (CCS[i]['pc'] > 0) and (F35[i]['pc'] < 0) and (F40[i]['pc'] < 0):
                A_B.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , F35[i]['pc'] , F40[i]['pc']))
            elif (CCS[i]['pc'] < 0) and (F35[i]['pc'] > 0) and (F40[i]['pc'] > 0):
                B_A.append((CCS[i]['pc'] , NT5[i]['pc'] , NT6[i]['pc'] , F35[i]['pc'] , F40[i]['pc']))
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
    ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'fuchsia','linewidth':2},
            medianprops={'color':'fuchsia','linewidth':2},
            capprops={'color':'fuchsia','linewidth':2},
            whiskerprops={'color':'fuchsia','linewidth':2})
    
    ax.plot([0.5,5.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 , 4 , 5])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F35' , 'F40' ] , fontsize = 28)
    ax.set_xlim((0.5 , 5.5))
    ax.set_ylim((-0.1 , 0.1))
    return fig
             
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
                   
                    
CCS = Load_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt')
NT5 = Load_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT5_Traditonal_PC_200K_Compartment_200K.txt')
NT6 = Load_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT6_Traditonal_PC_200K_Compartment_200K.txt')
F35 = Load_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F35_Traditonal_PC_200K_Compartment_200K.txt')
F40 = Load_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F40_Traditonal_PC_200K_Compartment_200K.txt')




A_B , B_A = PC_changed(CCS , NT5 , NT6 , F35 , F40)

data1 = [A_B['CCs'] , A_B['NT5'] , A_B['NT6'] , A_B['F35'] , A_B['F40']]
data2 = [B_A['CCs'] , B_A['NT5'] , B_A['NT6'] , B_A['F35'] , B_A['F40']]

fig1 = Box_plot(data1)
fig2 = Box_plot(data2)
run_Plot(fig1 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\Fig2\\compartment_A_B_200K.pdf')
run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\Fig2\\compartment_B_A_200K.pdf')
