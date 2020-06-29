# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 16:26:45 2019

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
                
 
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
    
    ax.spines['right'].set_lw(0.5)
    ax.spines['right'].set_color('#b3b3b3')
    ax.spines['right'].set_alpha(0.9)

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
    


def Get_Loop_Strength(inter,start,end,Len):
    """
    """
    A = inter[(start - 1 * Len) : (start + 2 * Len), (end - 1 * Len): (end + 2 * Len)]
    B = inter[(start - 4 * Len) : (start - 1 * Len), (end + 2 * Len): (end + 5 * Len)]
    n_A = len(np.nonzero(A)[0])
    n_B = len(np.nonzero(B)[0])
    if n_A != 0:
        A_S = A.sum() / n_A
    else:
        A_S = 0
    if n_B != 0:
        B_S = B.sum() / n_B
    else:
        B_S = 0.0001
    return A_S /B_S

def Get_loops(LoopSource):
    """
    cell : String ,The name of cells  e.g.'fESC' , 'ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loops = []
    Loop = np.loadtxt(LoopSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append(i)
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops


#--------------------------------------------------Files----------------------------------------------------------------    
        
R = 20000
res = '20K'

CCS_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_' + res + '/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_' + res + '/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_' + res + '/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_' + res + '/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}
    

OutFolder = '/public/home/xxli/data/BDF1_New/HiC/Plot/selected_interval_20K_loop_0.1_ww4'
LoopFolder = '/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.1_ww4'

chrom = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']


size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})

selected_interval = [('3' , 33760000 , 35680000 , 'Sox2') , ('4' , 54680000 , 56840000 , 'Klf4') , ('17' , 35280000 , 36360000 , 'Oct4') , 
                     ('6' , 121400000 , 123320000 , 'Nanog') , ('12' , 86080000 , 87120000 , 'Esrrb') , ('2' , 106900000 , 108000000 , 'Fig2')]


for c in cell:
    OutFil = 'Heatmap_' + c + '_selected_20K_rm2_4.pdf'
    HiCData = HiC_Data[c]
    LoopFil = 'Cluster_' + c + '_loops20K_rm2_4.txt'
    LoopFil_raw = c + '_loops_20K.txt'
    LoopSource = os.path.join(LoopFolder , LoopFil)
    LoopSource_raw = os.path.join(LoopFolder , LoopFil_raw)
    LoopData = Get_loops(LoopSource)
    LoopData_raw = Get_loops(LoopSource_raw)
    pp = PdfPages(os.path.join(OutFolder , OutFil))
    for i in selected_interval:
        lib = HiCData[i[0]]
        startHiC = i[1] // R
        endHiC = i[2] // R
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        if nonzero.shape[0] == 0:
            break
        vmax = np.percentile(nonzero, 95)
        loops = LoopData[LoopData['chr'] == i[0]]
        loops_raw = LoopData_raw[LoopData_raw['chr'] == i[0]]
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labels = [properU(p) for p in pos]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlabel('chr' + i[0] + '_' + i[-1], fontsize=15,labelpad=10)
        ## loop data
        mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
        extract = loops[mask]
        mask_raw = (loops_raw['start'] >= startHiC * R) & (loops_raw['end'] <= endHiC * R)
        extract_raw = loops_raw[mask_raw] 
        for p in extract:
            x = p['start'] // R - startHiC
            y = p['end'] // R - startHiC
            ax.scatter(y + 0.5 ,x + 0.5 , color = '', edgecolors = 'b', s = 100)
        for p in extract_raw:
                x = p['end'] // R - startHiC
                y = p['start'] // R - startHiC
                ax.scatter(y + 0.5,x + 0.5 , marker = 's' , color = 'b', edgecolors = '', s = 2)
#            strength = Get_Loop_Strength(matrix , y  , x , 1)
#            ax.text(y + 1.5 , x + 1 , str(round(strength,2)) , size = 15)
                
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        caxis_H(ax)
            
        ## Colorbar
        ax = fig.add_axes([Left + 0.5 , HB - 0.1 , 0.1 , 0.035])
        cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
        cbar.set_ticks([0 , int(vmax/2) , int(vmax)])

        pp.savefig(fig)
        plt.close(fig)
        ##Loopscatter
        
    pp.close()