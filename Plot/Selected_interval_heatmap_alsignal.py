# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:56:14 2019

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
import pyBigWig
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
    HiC_Data = np.load(HiC_fil)
    Lib = HiC_Data['Matrix'][()]
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
    
def bigwig_10bp_Plot(fil,chro,start,end,fig,location,color,label):
    """
    """
    sig_type = np.dtype({'names':['start' , 'end' , 'score'],
                      'formats':[np.int , np.int , np.float]})
    bw = pyBigWig.open(fil)
    bw = bw.intervals(chro, start, end)
    tmp_data = np.array(list(bw) , dtype = sig_type)
    bin_size = (end - start) // 10 + 1
    sig_data = np.zeros((bin_size,))
    for line in tmp_data:
        s = line['start'] // 10 - start // 10
        e = line['end'] // 10 - start // 10
        for i in range(s,e):
            if i >= 0 and i < bin_size: 
                sig_data[i] += line['score']
            else:
                pass
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
    ax.set_xlim((0 , len(sig_data)))
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)
    
    

    



        
#--------------------------------------------------Files----------------------------------------------------------------


OutFolder = '/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4'

cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3 , 'CH12':4}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC' , 'CH12']
R = 20000
res = '20K'

#selected_interval = [('3' , 33760000 , 35680000 , 'Sox2') , ('4' , 55200000 , 57000000 , 'Klf4') , 
#                     ('17' , 35280000 , 36360000 , 'Oct4') , ('6' , 121400000 , 123320000 , 'Nanog') , 
#                     ('12' , 86080000 , 87120000 , 'Esrrb') , ('2' , 106900000 , 108000000 , 'Fig2')]
selected_interval = [('4' , 55200000 , 57000000 , 'Klf4')]


size = (12, 12)
Left = 0.25 ; HB = 0.05 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    
loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})
    
    
f = open('/public/home/xxli/data/ref/haplotype/mm10.txt' , 'r')
mm ={}
for i in f:
    i = i.split()
    mm[i[0]] = int(i[1])
    

    



tmp = {}

# Data process
CCS_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CH12_DI = DI_Calling('/public/home/xxli/data/literature_data/cell_3D_map/workspace/Merge_6reps/CH12_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
DI_Data = {'CCS':CCS_DI,
           'fESC':fESC_DI,
           'NT5':NT5_DI,
           'NT6':NT6_DI,
           'CH12':CH12_DI}


CCS_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
CH12_Lib = np.load('/public/home/xxli/data/literature_data/cell_3D_map/workspace/Merge_6reps/CH12_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CH12_Lib = CH12_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'CH12':CH12_Lib}

CCS_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/CCS_compartment_200K.txt')
NT5_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT5_compartment_200K.txt')
NT6_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT6_compartment_200K.txt')
fESC_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/fESC_compartment_200K.txt')
CH12_PCA = pca_To_20K('/public/home/xxli/data/literature_data/cell_3D_map/workspace/Merge_6reps/CH12_200K/Paternal_Compartment/Paternal_Compartment_Compartment_200K.txt')
PCA_Data = {'CCS':CCS_PCA,
            'fESC':fESC_PCA,
            'NT5':NT5_PCA,
            'NT6':NT6_PCA,
            'CH12':CH12_PCA}

    
CCS_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_CCS_loops20K_3.txt')
NT5_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_NT5_loops20K_3.txt')
NT6_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_NT6_loops20K_3.txt')
fESC_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_fESC_loops20K_3.txt')
CH12_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_CH12_loops20K_6reps_3.txt')
Loop_Data = {'CCS':CCS_Loop,
            'NT5':NT5_Loop,
            'NT6':NT6_Loop,
            'fESC':NT5_Loop,
            'CH12':CH12_Loop}
            
CCS_RNA = '/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/CCS_RNA_10bp.bw'
NT5_RNA = '/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/NT5_RNA_10bp.bw'
NT6_RNA = '/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/NT6_RNA_10bp.bw'
fESC_RNA = '/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/fESC_RNA_10bp.bw'
CH12_RNA = '/public/home/xxli/data/literature_data/cell_3D_map/RNA/CH12_RNA.bw'
RNA_Data = {'CCS':CCS_RNA,
            'NT5':NT5_RNA,
            'NT6':NT6_RNA,
            'fESC':fESC_RNA,
            'CH12':CH12_RNA}


CCS_ATAC = '/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/CCS_ATAC.bw'
NT5_ATAC = '/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/NT5_ATAC.bw'
NT6_ATAC = '/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/NT6_ATAC.bw'
fESC_ATAC = '/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/fESC_ATAC.bw'
CH12_DNase = '/public/home/xxli/data/literature_data/cell_3D_map/DNase/CH12_DNase.bw'
ATAC_Data = {'CCS':CCS_ATAC,
            'NT5':NT5_ATAC,
            'NT6':NT6_ATAC,
            'fESC':fESC_ATAC,
            'CH12':CH12_DNase}

pp1 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4/Selected_20K_CCS_alsignal.pdf')
pp2 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4/Selected_20K_NT5_alsignal.pdf')
pp3 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4/Selected_20K_NT6_alsignal.pdf')
pp4 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4/Selected_20K_fESC_alsignal.pdf')
pp5 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4/Selected_20K_CH12_alsignal.pdf')

for i in selected_interval:
    print i
    g = i[0]
    startHiC = i[1] // R 
    endHiC = i[2] // R 
    for c in cell:
        lib = HiC_Data[c][g]
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
#        ax1.set_xlabel(i[-1] , fontsize = 30 )
        
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
        ax2 = fig.add_axes([Left + width + 0.015 , HB , 0.035 , 0.1])
        fig.colorbar(sc,cax = ax2, orientation='vertical' , ticks = [matrix.min() , vmax])

        ##gene Tracks
        ax3 = fig.add_axes([Left , HB + HH , width , 0.04])
        caxis_S_horizontal(ax3, 'black')
        
        if c == 'CH12':
            chro = 'chr' + g
        else:
            chro = g 
        ##RNA Tracks
        location = [Left , HB + HH + 0.05, width , 0.075]
        ax4 = fig.add_axes(location)
        bigwig_10bp_Plot(RNA_Data[c] , chro , i[1] , i[2] , fig , location , 'fuchsia' , 'RNA')
        ##ATAC Tracks
        location = [Left , HB + HH + 0.05 + 0.075, width , 0.075]
        ax5 = fig.add_axes(location)
        bigwig_10bp_Plot(ATAC_Data[c] , chro , i[1] , i[2] , fig , location , 'mediumblue' , 'ATAC')       
        ##DI Tracks
        ax6 = fig.add_axes([Left , HB + HH + 0.05 + 0.075 + 0.075, width , 0.075])
        DI_index , DI = UpdateDI(DI)
        ax6.fill_between(np.arange(DI.size) , DI , where = DI >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax6.fill_between(np.arange(DI.size) , DI , where = DI <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax6.set_ylabel('DI',fontsize = 15,rotation = 'horizontal',labelpad = 50)
        ax6.set_xlim((cxlim[0] , cxlim[1] - 1))
        caxis_S_horizontal(ax6, 'black')
        ##PCA Tracks
        PCA_index , PCA = UpdateDI(PCA)
        ax7 = fig.add_axes([Left, HB + width + 0.05 + 0.075 + 0.075 + 0.075, width , 0.075])
        ax7.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax7.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax7.set_xlim(0 , PCA_index.max())
        ax7.set_ylabel('PC1',fontsize = 15,rotation = 'horizontal',labelpad = 50)
        caxis_S_horizontal(ax7, 'black')
        
        tmp[c] = fig
    
    
    ccs_fig = tmp['CCS']
    fesc_fig = tmp['fESC']
    nt5_fig = tmp['NT5']
    nt6_fig = tmp['NT6']
    ch12_fig = tmp['CH12']
    
    figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig,ch12_fig]
    standard_axes_lim(3,*figs)
    standard_axes_lim(4,*figs)
    standard_axes_lim(5,*figs)
    standard_axes_lim(6,*figs)

    pp1.savefig(ccs_fig)
    pp2.savefig(nt5_fig)
    pp3.savefig(nt6_fig)
    pp4.savefig(fesc_fig)
    pp5.savefig(ch12_fig)
    plt.close(ccs_fig)
    plt.close(fesc_fig)
    plt.close(nt5_fig)
    plt.close(nt6_fig)
    plt.close(ch12_fig)

pp1.close()
pp2.close()
pp3.close()
pp4.close()
pp5.close()