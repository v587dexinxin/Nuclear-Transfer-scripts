# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:58:14 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Our Own Color Map
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')

#---------------------------------------------------Functions-------------------------------------------------------------
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
    
def caxis_S_vertical(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'y', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
def caxis_S_horizontal(ax):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
        
#-----------------------------------------------Files-------------------------------------------------------------------------    
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

cell = ['CCS' , 'MEF' , 'NT5' , 'NT6' , 'fESC' , 'IPS_P3' , 'E14']
cells = {'CCS':0 , 'MEF':1 , 'NT5':2 , 'NT6':3 , 'fESC':4 , 'IPS_P3':5 , 'E14':6}
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
R = 200000
res = '200K'

Outfolder = 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\'

CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = sig_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = sig_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = sig_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = sig_type)
fESC = fESC[fESC['chr'] != 'X']


MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = sig_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = sig_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = sig_type)
E14 = E14[E14['chr'] != 'X']


PC_all = {'CCS':CCS['score'] , 'MEF':MEF['score'] , 'NT5':NT5['score'] , 'NT6':NT6['score'] , 
          'fESC':fESC['score'] , 'IPS_P3':IPS_P3['score'] , 'E14':E14['score']}
    


cor_matrix = np.zeros((7,7))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(PC_all[i] , PC_all[j])[0][1] , 3)
        
pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\Compartment_Correlation_heatmap_200K.pdf')
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
size_axes = [Left, HB, width, HH]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
vmax = 0.99
vmin = 0.83
im = ax.imshow(cor_matrix , vmax=vmax , vmin = vmin , cmap = my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of Compartment between different cells',fontsize=25)

for i in range(len(x)):
    ax.plot([-0.5 , 6.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 6.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 20)

##Colorbar
ax = fig.add_axes([Left + width + 0.035 , HB , 0.035 , 0.1])
cbar = fig.colorbar(im,cax = ax, orientation='vertical')
cbar.set_ticks([vmin , vmax])

pp.savefig(fig)
pp.close()

##hierarchy_cluster
data = np.vstack((PC_all['CCS'] , PC_all['MEF'] , PC_all['NT5'] , PC_all['NT6'] , PC_all['fESC'], PC_all['IPS_P3'] , PC_all['E14']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)