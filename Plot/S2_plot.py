# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 15:30:44 2019

@author: Administrator
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
    
    


    
def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'x', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.tick_params(axis = 'y', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0.9)
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
        s_f.axes[N].set_ylim(np.round(min(lim_1), 2 ),np.round(max(lim) , 2) + 0.01)
        s_f.axes[N].set_yticks([np.round(min(lim_1) * 0.5 , 2) , np.round(max(lim) * 0.5 , 2)])
        
        
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

def Sig_To_1K(fil):
    """
    """
    sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    ATACData = np.loadtxt(fil , dtype=sig_type)
    
    chroms = set(ATACData['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = ATACData[ATACData['chr'] == c]
        max_ = tmp_data['end'].max()
        bin_size = max_ // 1000 + 1
        New_Data[c] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 1000
            end = line['end'] // 1000
            for i in range(start,end):
                New_Data[c][i] += line['score']
    
    return New_Data


    
    
def Sig_Plot(data,start,end,chro,fig,location,color,label):
    """
    """
    tmp = data[chro]
    start = start * 20
    end = end * 20
    sig_data = tmp[start:end]
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
#    ax.set_ylabel(label,fontsize = 40,rotation = 'horizontal',labelpad = 50)
    ax.set_xlim([0 , len(sig_data)])
    ax.set_ylim([sig_data.min() , sig_data.max()])
    caxis_S(ax,'black')
    
        
#--------------------------------------------------Files----------------------------------------------------------------


OutFolder = '/public/home/xxli/data/BDF1_New/HiC/Plot/Fig4'

cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3 , 'CH12':4}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC' ]
R = 20000
res = '20K'

selected_interval = [('3' , 33760000 , 35680000 , 'Sox2') , ('4' , 55200000 , 57000000 , 'Klf4') , 
                     ('17' , 35280000 , 36360000 , 'Oct4') , ('6' , 122300000 , 123320000 , 'Nanog') , 
                     ('12' , 86080000 , 87120000 , 'Esrrb') , ('2' , 106900000 , 108000000 , 'Fig2')]



size = (12, 12)
Left = 0.25 ; HB = 0.1 ; width = 0.5 ; HH = 0.5

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
    

    





# DI Data process
CCS_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_DI = DI_Calling('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
DI_Data = {'CCS':CCS_DI,
           'fESC':fESC_DI,
           'NT5':NT5_DI,
           'NT6':NT6_DI}


# TAD Data process
CCS_Tad = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/HiTAD/respective/stable_600K/CCS_40K_allreps_union_small_domain.txt' , dtype = tad_type)
fESC_Tad = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/HiTAD/respective/stable_600K/CCS_40K_allreps_union_small_domain.txt' , dtype = tad_type)
NT5_Tad = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/HiTAD/respective/stable_600K/CCS_40K_allreps_union_small_domain.txt' , dtype = tad_type)
NT6_Tad = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/HiTAD/respective/stable_600K/CCS_40K_allreps_union_small_domain.txt' , dtype = tad_type)
Tad_Data = {'CCS':CCS_Tad,
           'fESC':fESC_Tad,
           'NT5':NT5_Tad,
           'NT6':NT6_Tad}

# HiC Data process
CCS_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/CCS_workspace/CCS_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/fESC_workspace/fESC_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT5_workspace/NT5_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/NT6_workspace/NT6_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}

# PCA DATA process
CCS_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/CCS_compartment_200K.txt')
NT5_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT5_compartment_200K.txt')
NT6_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT6_compartment_200K.txt')
fESC_PCA = pca_To_20K('/public/home/xxli/data/BDF1_New/HiC/Compartment/fESC_compartment_200K.txt')
PCA_Data = {'CCS':CCS_PCA,
            'fESC':fESC_PCA,
            'NT5':NT5_PCA,
            'NT6':NT6_PCA}

#Loop Data process
CCS_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_CCS_loops20K_3.txt')
NT5_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_NT5_loops20K_3.txt')
NT6_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_NT6_loops20K_3.txt')
fESC_Loop = Get_loops('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_fESC_loops20K_3.txt')
Loop_Data = {'CCS':CCS_Loop,
            'NT5':NT5_Loop,
            'NT6':NT6_Loop,
            'fESC':NT5_Loop}

# ATAC data process
CCS_ATAC = Sig_To_1K('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/CCS_ATAC_1K.bedgraph')
fESC_ATAC = Sig_To_1K('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/fESC_ATAC_1K.bedgraph')
NT5_ATAC = Sig_To_1K('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT5_ATAC_1K.bedgraph')
NT6_ATAC = Sig_To_1K('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT6_ATAC_1K.bedgraph')
ATAC_Data = {'CCS':CCS_ATAC,
             'fESC':fESC_ATAC,
             'NT5':NT5_ATAC,
             'NT6':NT6_ATAC}


#RNA Data process
CCS_RNA = Sig_To_1K('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bedgraph_1K/CCS_RNA_1K.bedgraph')
fESC_RNA = Sig_To_1K('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bedgraph_1K/fESC_RNA_1K.bedgraph')
NT5_RNA = Sig_To_1K('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bedgraph_1K/NT5_RNA_1K.bedgraph')
NT6_RNA = Sig_To_1K('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bedgraph_1K/NT6_RNA_1K.bedgraph')
RNA_Data = {'CCS':CCS_RNA,
            'fESC':fESC_RNA,
            'NT5':NT5_RNA,
            'NT6':NT6_RNA}


tmp = {}
pp1 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/S2/Selected_20K_CCS_loops.pdf')
pp2 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/S2/Selected_20K_NT5_loops.pdf')
pp3 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/S2/Selected_20K_NT6_loops.pdf')
pp4 = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Plot/S2/Selected_20K_fESC_loops.pdf')

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
        tads = Tad_Data[c][Tad_Data[c]['chr'] == g]
        
        
        #=============HeatMap + colorbar + DI=================================
        
        fig = plt.figure(figsize = size)
        ax1 = fig.add_axes([Left  , HB , width , HH])
        sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax1.get_xlim()
        cylim = ax1.get_ylim()
        ## Ticks and Labels
        ticks = [0 , len(matrix)]
#        pos = [((startHiC + t) * R) for t in ticks]
#        labels = [str(np.round(startHiC * R /1000000 , 1)) + 'M' , str(np.round(endHiC * R /1000000 , 1)) + 'M']
#        ax1.set_xticks(ticks)
#        ax1.set_xticklabels(labels)
#        ax1.set_yticks(ticks)
#        ax1.set_yticklabels(labels, rotation = 'vertical')
#        ax1.set_xlabel(i[-1] , fontsize = 30 )
        
#        ## loop data
        mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
        extract = loops[mask]
        for p in extract:
            x = p['start'] // R - startHiC
            y = p['end'] // R - startHiC
            ax1.scatter(y + 0.5 ,x + 0.5 , color = '', edgecolors = 'black', s = 200)
        ## Domain Boundaries
        mask = (tads['end'] > startHiC * R ) & (tads['start'] < endHiC * R )
        extract = tads[mask]
        pairs = [(bi['start']//R - startHiC, bi['end']//R - startHiC) for bi in extract]
        for corner in pairs:
            if (corner[0] <= 0):
                ax1.plot([0, corner[1]] , [corner[1] , corner[1]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
                ax1.plot([corner[1] , corner[1]] , [0, corner[1]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
            elif (corner[1] >= len(matrix)):
                ax1.plot([corner[0] , corner[0]] , [corner[0] , len(matrix)] , color = 'black', linewidth = 2, linestyle = (0, (5, 10)))
                ax1.plot([corner[0] , len(matrix)] , [corner[0] , corner[0]] , color = 'black', linewidth = 2, linestyle = (0, (5, 10)))
            else:
                ax1.plot([corner[0] , corner[0]] , [corner[0] , corner[1]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
                ax1.plot([corner[0] , corner[1]] , [corner[0] , corner[0]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
                ax1.plot([corner[0] , corner[1]] , [corner[1] , corner[1]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
                ax1.plot([corner[1] , corner[1]] , [corner[0] , corner[1]] , color = 'black' , linewidth = 2, linestyle = (0, (5, 10)))
            
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)                    
        caxis_H(ax1)
        ## Colorbar
        ax2 = fig.add_axes([Left + 0.2 , HB - 0.06 , 0.1 , 0.035])
        fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0 , int(vmax)])
        caxis_colorbar(ax2)

        ##DI Tracks
        ax3 = fig.add_axes([Left , HB + HH , width , 0.075])
        DI_index , DI = UpdateDI(DI)
        ax3.fill_between(np.arange(DI.size) , DI , where = DI >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax3.fill_between(np.arange(DI.size) , DI , where = DI <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
#        ytick = [DI.min() + 1 , DI.max() - 1]
#        ax3.set_yticks(ytick)
#        ax3.set_ylabel('DI',fontsize=20,)
        ax3.set_xlim((cxlim[0] , cxlim[1] - 1))
        caxis_S(ax3, 'black')
        ##PCA Tracks
        PCA_index , PCA = UpdateDI(PCA)
        ax4 = fig.add_axes([Left, HB + HH + 0.075, width , 0.075])
        ax4.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax4.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax4.set_xlim(0 , PCA_index.max())
#        ytick = [round(sig.min() * 0.6667 , 2) , 0.00 , round(sig.max() * 0.6667, 2)]
#        ax.set_yticks(ytick)
#        ax4.set_ylim((PCA.min() , PCA.max()))
        caxis_S(ax4, 'black')

        #ATAC
        location1 = [Left, HB + HH + 0.15 , width , 0.075]
        color = 'blue'
        Sig_Plot(data = ATAC_Data[c],
                 start = startHiC,
                 end = endHiC,
                 chro = g,
                 fig = fig,
                 location = location1,
                 color = color,
                 label = 'ATAC')
        #RNA
        location2 = [Left , HB + HH + 0.225, width , 0.075]
        color = 'fuchsia'
        Sig_Plot(data = RNA_Data[c],
                 start = startHiC,
                 end = endHiC,
                 chro = g,
                 fig = fig,
                 location = location2,
                 color = color,
                 label = 'RNA')
        
        
        
            
            
        tmp[c] = fig
    
    
    ccs_fig = tmp['CCS']
    fesc_fig = tmp['fESC']
    nt5_fig = tmp['NT5']
    nt6_fig = tmp['NT6']

    
    figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig]
    standard_axes_lim(2,*figs)
    standard_axes_lim(3,*figs)
    standard_axes_lim(4,*figs)
    standard_axes_lim(5,*figs)

    pp1.savefig(ccs_fig)
    pp2.savefig(nt5_fig)
    pp3.savefig(nt6_fig)
    pp4.savefig(fesc_fig)
    plt.close(ccs_fig)
    plt.close(fesc_fig)
    plt.close(nt5_fig)
    plt.close(nt6_fig)


pp1.close()
pp2.close()
pp3.close()
pp4.close()
