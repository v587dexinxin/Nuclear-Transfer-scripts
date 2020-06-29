# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 21:15:52 2019

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

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

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
    
def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

def add_ax(sig , loc , cell):
    PCA_index1 , PCA1 = UpdateDI(sig)
    ax = fig.add_axes(loc)
    ax.fill_between(PCA_index1 , PCA1 , where = PCA1 >= 0 , facecolor = 'gold' , edgecolor = 'none' )
    ax.fill_between(PCA_index1 , PCA1 , where = PCA1 <= 0 , facecolor = 'midnightblue' , edgecolor = 'none' )
    ytick = [round(sig.min() * 0.6667 , 2) , 0.00 , round(sig.max() * 0.6667, 2)]
    ax.set_yticks(ytick)
    ax.set_ylim((sig.min() * 1.1 , sig.max()* 1.1))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)
    return ax
        
#-----------------------------------------------Files-------------------------------------------------------------------------    

chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
R = 200000
res = '200K'

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

PCFolder = 'H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new'
Outfolder = 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\'


MEF = np.loadtxt(os.path.join(PCFolder , 'MEF_compartment_200K.txt') , dtype = sig_type)
IPS_P3 = np.loadtxt(os.path.join(PCFolder , 'IPS_P3_compartment_200K.txt') , dtype = sig_type)
E14 = np.loadtxt(os.path.join(PCFolder , 'E14_compartment_200K.txt') , dtype = sig_type)



#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 10)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7


OutFil = 'IPSC_Compartment_200K.pdf'
pp = PdfPages(os.path.join(Outfolder , OutFil))


for g in chrom:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]

    fig = plt.figure(figsize = size)
    sig1 = tmp_MEF['score']
    sig2 = tmp_IPS_P3['score']
    sig3 = tmp_E14['score']
    ##PCA Tracks
    ax = add_ax(sig1 , [Left, HB , width , 0.1] , 'MEF')
    xtick = list(np.linspace(0 , len(sig1) , 5).astype(float))
    pos = [(t * R) for t in xtick]
    labels = [properU(p) for p in pos]
    ax.set_xticks(xtick)
    ax.set_xticklabels(labels)
    ax.set_xlabel(g , fontsize=20)
    ax = add_ax(sig2 , [Left, HB + 0.1 , width , 0.1] , 'IPS_P3')
    caxis_S_horizontal(ax)
    ax = add_ax(sig3 , [Left, HB + 0.2, width , 0.1] , 'E14')
    caxis_S_horizontal(ax)

        
        
        
        
    pp.savefig(fig)
    plt.close(fig)
pp.close()
    


