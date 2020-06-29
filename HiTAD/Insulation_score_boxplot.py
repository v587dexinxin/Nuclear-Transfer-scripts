# -*- coding: utf-8 -*-
"""
Created on Sun Nov 03 21:19:31 2019

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

IS_type = np.dtype({'names':['chr' , 'IS'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float]})
 

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
#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F40' ] , fontsize = 28)
    ax.set_xlim((0.5 , 4.5))
#    ax.set_ylim((-0.1 , 0.1))
    return fig
             
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
                   
                    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\Insulation_score\\CCS_Insulation_score_chrX_40K.txt' , dtype = IS_type , skiprows = 1 , usecols = (0 , 2))
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\Insulation_score\\NT5_Insulation_score_chrX_40K.txt' , dtype = IS_type , skiprows = 1 , usecols = (0 , 2))
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\Insulation_score\\NT6_Insulation_score_chrX_40K.txt' , dtype = IS_type , skiprows = 1 , usecols = (0 , 2))
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\Insulation_score\\fESC_Insulation_score_chrX_40K.txt' , dtype = IS_type , skiprows = 1 , usecols = (0 , 2))

CCS = CCS[CCS['chr'] == 'X']
NT5 = NT5[NT5['chr'] == 'X']
NT6 = NT6[NT6['chr'] == 'X']
fESC = fESC[fESC['chr'] == 'X']

#chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

chroms = ['X']


data = [CCS['IS'] , NT5['IS'] , NT6['IS'] , fESC['IS']]


fig = Box_plot(data)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\ChrX_Insulation_score_boxplot.pdf')
