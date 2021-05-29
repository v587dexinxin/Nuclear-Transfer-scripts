# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 15:59:04 2020

@author: xxli
"""



from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['lightcyan' , '#FFFFFF' , 'lavender'])
my_cmap.set_bad('#D3D3D3')
                
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 25 , pad = 5 )
                
cell = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
# chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
# chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18}
PCFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new'
#OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\'
size = (12, 12)   
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U4',np.float]})

pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\correlation\\Compartment_Correlation_replicates_scatter_200K_1.pdf')
for c in cell:
    PC_R1 = np.loadtxt(os.path.join(PCFolder , c + '_R1_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
    PC_R1 = PC_R1[PC_R1['chr'] != 'X']
    PC_R2 = np.loadtxt(os.path.join(PCFolder , c + '_R2_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
    PC_R2 = PC_R2[PC_R2['chr'] != 'X']
    cor = round(np.corrcoef(PC_R1['score'] , PC_R2['score'])[0][1],5)
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(PC_R1['score'] , PC_R2['score'] , alpha = 0.8 , c = 'red')
    ax.set_xlim(-0.13,0.13)
    ax.set_ylim(-0.13,0.13)
    ax.set_xticks([-0.13 , 0.13])
    ax.set_xticklabels([-0.13 , 0.13] , size = 25 )
    ax.set_yticks([-0.13 , 0.13])
    ax.set_yticklabels([-0.13 , 0.13] , size = 25 )
    ax.text(-0.1 , 0.1 , 'PCC = ' + str(np.round(cor,3)) , size = 40 )
    ax.set_xlabel(c + '_R1' , size = 50 , labelpad = 0.001)
    ax.set_ylabel(c + '_R2' , size = 50 , labelpad = 0.001)
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
#    caxis_H(ax)
    pp.savefig(fig)
pp.close()