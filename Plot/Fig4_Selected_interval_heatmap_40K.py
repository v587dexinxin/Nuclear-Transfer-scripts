# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:09:05 2020

@author: han-luo
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
                

def binbias(up, down):
    bias = 0
    zeromask = (up != 0) & (down != 0)
    if  zeromask.sum() <= 1:
        return bias
        
    upmean = up.mean(); downmean = down.mean()
    SD_1 = np.sum((up - upmean) ** 2) / (up.size * (up.size - 1))
    SD_2 = np.sum((down - downmean) ** 2) / (down.size * (down.size - 1))
    SD_pool = np.sqrt(SD_1 + SD_2)
    if SD_pool != 0:
        bias = (upmean - downmean) / SD_pool

    return bias
    
    
def CalDI(matrix):
    n = 0
    DIs = []
    shape = matrix.shape
    for j in matrix[:shape[0]]:
        if n <20:
            DIs.append(0)
        elif n >= 20 and n <= int(shape[0])-21:
            if len(j[j!=0])/int(shape[0]) < 0.05:
                bias = 0
            else:
                up = j[n-20:n] 
                down = j[n+1:n+21]
                bias = binbias(up,down) 
            DIs.append(bias)
        else:
            DIs.append(0)
        n += 1
    return DIs


def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 , labelsize = 23)
    
    
    
def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
    

def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

    
def caxis_DI(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
    

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
    
def DI_Calling(HiC_fil):
    """
    """
    Lib = np.load(HiC_fil , allow_pickle=True)
    Lib = Lib['Matrix'][()]
    DI_Data = {}
    for chro in Lib.keys():
        Matrix = Lib[chro]
        DI_score = CalDI(Matrix)
        DI_Data[chro] = DI_score
    
    return DI_Data



        
def pca_To_20K(fil):
    """
    """
    pca_type = np.dtype({'names':['chr' , 'PCA'],
                     'formats':['S4' , np.float]})
    PCA_Data = np.loadtxt(fil , dtype=pca_type)
    
    chroms = set(PCA_Data['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = PCA_Data[PCA_Data['chr'] == c]
        New_Data[c] = []
        for i in tmp_data:
            New_Data[c].extend([i['PCA']] * 10)
            
    
    return New_Data


def pca_To_40K(fil):
    """
    """
    pca_type = np.dtype({'names':['chr' , 'PCA'],
                     'formats':['S4' , np.float]})
    PCA_Data = np.loadtxt(fil , dtype=pca_type)
    
    chroms = set(PCA_Data['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = PCA_Data[PCA_Data['chr'] == c]
        New_Data[c] = []
        for i in tmp_data:
            New_Data[c].extend([i['PCA']] * 5)
            
    
    return New_Data
    
    


def standard_axes_lim(N,*figs):
    """
    """    
    lim = []
    lim_1 = []
    for s_f in figs:
        lim.append(s_f.axes[N].get_ylim()[1])
        lim_1.append(s_f.axes[N].get_ylim()[0])
    
    for s_f in figs:
        s_f.axes[N].set_ylim(np.round(min(lim_1), 2 ) - 0.01,np.round(max(lim) , 2) + 0.01)
        
        
        
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



cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
R = 40000
res = '40K'

selected_interval = [('4' , 55200000 , 57000000 , 'Klf4')]



size = (12, 12)
Left = 0.2 ; HB = 0.15 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    
loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})
    
    
f = open('G:\\data\\genome\\mm10.txt' , 'r')
mm ={}
for i in f:
    i = i.split()
    mm[i[0]] = int(i[1])
    

    



tmp = {}

# DI Data process
CCS_DI = DI_Calling('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_DI = DI_Calling('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT5_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_DI = DI_Calling('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT6_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_DI = DI_Calling('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\fESC_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')

DI_Data = {'CCS':CCS_DI,
           'NT5':NT5_DI,
           'NT6':NT6_DI,
           'fESC':fESC_DI}


CCS_Lib = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT5_Lib = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT5_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT6_Lib = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\NT6_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
fESC_Lib = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\fESC_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
CCS_Lib = CCS_Lib['Matrix'][()]
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = NT6_Lib['Matrix'][()]
fESC_Lib = fESC_Lib['Matrix'][()]

CCS_Lib = CCS_Lib['4']
NT5_Lib = NT5_Lib['4']
NT6_Lib = NT6_Lib['4']
fESC_Lib = fESC_Lib['4']

HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'fESC':fESC_Lib}

CCS_PCA = pca_To_40K('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt')
NT5_PCA = pca_To_40K('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt')
NT6_PCA = pca_To_40K('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt')
fESC_PCA = pca_To_40K('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt')

PCA_Data = {'CCS':CCS_PCA,
            'NT5':NT5_PCA,
            'NT6':NT6_PCA,
            'fESC':fESC_PCA}

    
CCS_Loop = Get_loops('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\Cluster_CCS_loops20K_3.txt')
NT5_Loop = Get_loops('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\Cluster_NT5_loops20K_3.txt')
NT6_Loop = Get_loops('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\Cluster_NT6_loops20K_3.txt')
fESC_Loop = Get_loops('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\Cluster_fESC_loops20K_3.txt')

Loop_Data = {'CCS':CCS_Loop,
            'NT5':NT5_Loop,
            'NT6':NT6_Loop,
            'fESC':fESC_Loop}

pp1 = PdfPages('H:\\Workspace_New\\data\\IPSC\\HiC\\Plot\\CCS_Klf4_40K.pdf')
pp2 = PdfPages('H:\\Workspace_New\\data\\IPSC\\HiC\\Plot\\NT5_Klf4_40K.pdf')
pp3 = PdfPages('H:\\Workspace_New\\data\\IPSC\\HiC\\Plot\\NT6_Klf4_40K.pdf')
pp4 = PdfPages('H:\\Workspace_New\\data\\IPSC\\HiC\\Plot\\fESC_Klf4_40K.pdf')


for i in selected_interval:
    print i
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R 
    for c in cell:
        lib = HiC_Data[c]
        DIData = DI_Data[c][g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        DI = np.array(DIData[startHiC : endHiC])
        PCAData = PCA_Data[c][g]
        PCA = np.array(PCAData[startHiC:endHiC])
        loops = Loop_Data[c][Loop_Data[c]['chr'] == g]
        
        
        #=============HeatMap + colorbar + DI=================================
        
        fig = plt.figure(figsize = size)
        ax1 = fig.add_axes([Left  , HB , width , HH])
        sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax1.get_xlim()
        cylim = ax1.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labels = [properU(p) for p in pos[:4]]
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels, rotation = 'horizontal')
        ax1.set_xlabel(c , fontsize = 30 )
        
        ## loop data
        mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
        extract = loops[mask]
        for p in extract:
            x = p['start'] // R - startHiC
            y = p['end'] // R - startHiC
            ax1.scatter(y + 0.5 ,x + 0.5 , color = '', edgecolors = 'black', s = 200)
            
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)                    
        caxis_H(ax1)
        ## Colorbar
#        ax2 = fig.add_axes([Left + 0.25 , HB - 0.08 , 0.1 , 0.035])
#        fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0 , int(vmax)])
#        caxis_colorbar(ax2)

        ##DI Tracks
        ax3 = fig.add_axes([Left , HB + HH , width , 0.1])
        DI_index , DI = UpdateDI(DI)
        ax3.fill_between(np.arange(DI.size) , DI , where = DI >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax3.fill_between(np.arange(DI.size) , DI , where = DI <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
#        ytick = [DI.min() + 1 , DI.max() - 1]
#        ax3.set_yticks(ytick)
#        ax3.set_ylabel('DI',fontsize=20,)
        ax3.set_xlim((cxlim[0] , cxlim[1] - 1))
        caxis_S_horizontal(ax3, 'black')
        ##PCA Tracks
        PCA_index , PCA = UpdateDI(PCA)
        ax4 = fig.add_axes([Left, HB + width + 0.1, width , 0.1])
        ax4.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax4.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax4.set_xlim(0 , PCA_index.max())
#        ytick = [round(sig.min() * 0.6667 , 2) , 0.00 , round(sig.max() * 0.6667, 2)]
#        ax.set_yticks(ytick)
#        ax4.set_ylim((PCA.min() , PCA.max()))
        caxis_S_horizontal(ax4, 'black')
        
        
            
            
        tmp[c] = fig
    
    
    ccs_fig = tmp['CCS']
    nt5_fig = tmp['NT5']
    nt6_fig = tmp['NT6']
    fesc_fig = tmp['fESC']
    
    
    figs = [ccs_fig,nt5_fig,nt6_fig,fesc_fig]
    standard_axes_lim(1,*figs)
    standard_axes_lim(2,*figs)

    pp1.savefig(ccs_fig)
    pp2.savefig(nt5_fig)
    pp3.savefig(nt6_fig)
    pp4.savefig(fesc_fig)
    
    plt.close(ccs_fig)
    plt.close(nt5_fig)
    plt.close(nt6_fig)
    plt.close(fesc_fig)


pp1.close()
pp2.close()
pp3.close()
pp4.close()