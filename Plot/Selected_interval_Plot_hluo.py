# -*- coding: utf-8 -*-
"""
Created on Sat Dec 01 10:43:06 2018

@author: xxli
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os


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
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

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
                   labelleft = True, labelright = False)
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

def Sig_Plot(data,start,end,chro,fig,location,color,label):
    """
    """
    tmp = data[chro]
    start = start * 40
    end = end * 40
    sig_data = tmp[start:end]
    ax = fig.add_axes(location)
    ax.fill_between(np.arange(len(sig_data)),sig_data, facecolor = color, edgecolor = 'none')
    ax.set_ylabel(label,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S(ax,color)


def standard_axes_lim(N,*figs):
    """
    """    
    lim = []
    for s_f in figs:
        lim.append(s_f.axes[N].get_ylim()[1])
    
    for s_f in figs:
        s_f.axes[N].set_ylim(0,max(lim))
        
#--------------------------------------------------Files----------------------------------------------------------------


TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_800K'
ATACFolder = 'D:\\Workspace_New\\data\\ATAC\\signal\\normalization\\bedgraph_1K'
RNAFolder = 'D:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K'
OutFolder = 'D:\\Workspace_New\\Plot\\Figure2\\somatic_gene'

cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
R = 40000
res = '40K'
interval = [('17' ,  35506018 , 35510772 , 'Oct4') , ('3' , 34650405 , 34652461 , 'Sox2') , ('4' , 55527143 , 55532466 , 'Klf4') , 
            ('15' , 61985391 , 61990374 , 'cMyc') , ('6' , 122707489 , 122714633 , 'Nanog') , ('12' , 86361117 , 86521628 , 'Esrrb') , 
            ('16' , 92601466 , 92826149 , 'Runx1') , ('19' , 5447703 , 5455945 , 'Fra1') , ('7' , 35119293 , 35121928 , 'Cebpa') , 
            ('2' , 167688915 , 167690418 , 'Cebpb') , ('12' ,  85473890 , 85477273 , 'Fos') , ('4' , 95049034 , 95052222 , 'Jun') , 
            ('9' , 32636221 , 32757820 , 'Ets1') , ('16' , 95702075 , 95721051 , 'Ets2') , ('7' , 112679318 , 112906807 , 'Tead1') , 
            ('7' , 45215753 , 45233644 , 'Tead2') , ('17' , 28331671 , 28350805 , 'Tead3') , ('6' , 128224288 , 128300823 , 'Tead4')]

size = (12, 12)
Left = 0.1 ; HB = 0.15 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    
f = open('G:\\data\\genome\\mm10.txt' , 'r')
mm ={}
for i in f:
    i = i.split()
    mm[i[0]] = int(i[1])
    

    



tmp = {}

# ATAC data process
CCS_ATAC = Sig_To_1K('D:/Workspace_New/data/ATAC/signal/normalization/bedgraph_1K/CCS_ATAC_1K.bedgraph')
fESC_ATAC = Sig_To_1K('D:/Workspace_New/data/ATAC/signal/normalization/bedgraph_1K/fESC_ATAC_1K.bedgraph')
NT5_ATAC = Sig_To_1K('D:/Workspace_New/data/ATAC/signal/normalization/bedgraph_1K/NT5_ATAC_1K.bedgraph')
NT6_ATAC = Sig_To_1K('D:/Workspace_New/data/ATAC/signal/normalization/bedgraph_1K/NT6_ATAC_1K.bedgraph')
ATAC_Data = {'CCS':CCS_ATAC,
             'fESC':fESC_ATAC,
             'NT5':NT5_ATAC,
             'NT6':NT6_ATAC}


#RNA Data process
CCS_RNA = Sig_To_1K('D:/Workspace_New/data/RNA/signal/normalization/bedgraph_1K/CCS_RNA_1K.bedgraph')
fESC_RNA = Sig_To_1K('D:/Workspace_New/data/RNA/signal/normalization/bedgraph_1K/fESC_RNA_1K.bedgraph')
NT5_RNA = Sig_To_1K('D:/Workspace_New/data/RNA/signal/normalization/bedgraph_1K/NT5_RNA_1K.bedgraph')
NT6_RNA = Sig_To_1K('D:/Workspace_New/data/RNA/signal/normalization/bedgraph_1K/NT6_RNA_1K.bedgraph')
RNA_Data = {'CCS':CCS_RNA,
            'fESC':fESC_RNA,
            'NT5':NT5_RNA,
            'NT6':NT6_RNA}

# DI Data process
CCS_DI = DI_Calling('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/CCS_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_DI = DI_Calling('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/fESC_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_DI = DI_Calling('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT5_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_DI = DI_Calling('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT6_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
DI_Data = {'CCS':CCS_DI,
           'fESC':fESC_DI,
           'NT5':NT5_DI,
           'NT6':NT6_DI}


CCS_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/CCS_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/fESC_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT5_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT6_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}

CCS_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/CCS_compartment_200K.txt')
NT5_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/NT5_compartment_200K.txt')
NT6_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/NT6_compartment_200K.txt')
fESC_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/fESC_compartment_200K.txt')
PCA_Data = {'CCS':CCS_PCA,
            'fESC':fESC_PCA,
            'NT5':NT5_PCA,
            'NT6':NT6_PCA}

    

pp1 = PdfPages('D:\\Workspace_New\\Plot\\Figure2\\Heatmap_respective_domain\\Selected_CCS_40K_Heatmap_stable_window_wise_800K_PCA.pdf')
pp2 = PdfPages('D:\\Workspace_New\\Plot\\Figure2\\Heatmap_respective_domain\\Selected_NT5_40K_Heatmap_stable_window_wise_800K_PCA.pdf')
pp3 = PdfPages('D:\\Workspace_New\\Plot\\Figure2\\Heatmap_respective_domain\\Selected_NT6_40K_Heatmap_stable_window_wise_800K_PCA.pdf')
pp4 = PdfPages('D:\\Workspace_New\\Plot\\Figure2\\Heatmap_respective_domain\\Selected_fESC_40K_Heatmap_stable_window_wise_800K_PCA.pdf')

for i in interval:
    print i
    g = i[0]
    startHiC = i[1] // R -25
    endHiC = i[2] // R + 25
    for c in cell:
        TadFil = c + '_40K_allreps.txt'
        TadSource = os.path.join(TadFolder , TadFil)
        TadData = np.loadtxt(TadSource , usecols = (0 , 1 , 2) ,dtype = tad_type)
        lib = HiC_Data[c][g]
        DIData = DI_Data[c][g]
        tads = TadData[TadData['chr'] == g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        DI = np.array(DIData[startHiC : endHiC])
        PCAData = PCA_Data[c][g]
        PCA = np.array(PCAData[startHiC:endHiC])
        
        
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
        ax1.set_xlabel(i[-1] , fontsize = 30 )
        
        ## Domain Boundaries
        mask = (tads['end'] > startHiC * R ) & (tads['start'] < endHiC * R )
        extract = tads[mask]
        new = []
        for x in extract:
            mask1 = (x['start'] <= extract['start'] ) & (x['end'] >= extract['end'])
            overlap1 = extract[mask1]
            if overlap1.size == 1:
                new.append(x)
            else:
                continue
        extract = new
        pairs = [(bi['start']//R - startHiC, bi['end']//R - startHiC) for bi in extract]
        for corner in pairs:
            if (corner[0] <= 0):
                ax1.plot([0, corner[1]], [corner[1], corner[1]], color = '#000000' , linewidth = 1)
                ax1.plot([corner[1], corner[1]], [0, corner[1]], color = '#000000' , linewidth = 1)
            elif (corner[1] >= interval):
                ax1.plot([corner[0], corner[0]], [corner[0], interval] , color = '#000000' , linewidth = 1)
                ax1.plot([corner[0], interval], [corner[0], corner[0]] , color = '#000000' , linewidth = 1)
            else:
                ax1.plot([corner[0], corner[0]], [corner[0], corner[1]] , color = '#000000' , linewidth = 1)
                ax1.plot([corner[0], corner[1]], [corner[0], corner[0]] , color = '#000000' , linewidth = 1)
                ax1.plot([corner[0], corner[1]], [corner[1], corner[1]] , color = '#000000' , linewidth = 1)
                ax1.plot([corner[1], corner[1]], [corner[0], corner[1]] , color = '#000000' , linewidth = 1)
                        
#b3b3b3
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)                    
        caxis_H(ax1)
        ## Colorbar
        ax2 = fig.add_axes([Left + 0.5 , HB - 0.08 , 0.1 , 0.035])
        fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0 , int(vmax/2) , int(vmax)])


        ##DI Tracks
        ax3 = fig.add_axes([Left + width , HB , 0.1 , HH])
        ax3.fill_betweenx(np.arange(DI.size) , DI , facecolor = 'blue' , edgecolor = 'none' )
        xtick = [DI.min() + 1 , DI.max() - 1]
        ax3.set_xticks(xtick)
        ax3.set_xlabel('DI',fontsize=20,)
        ax3.set_xlim((DI.min() , DI.max()))
        ax3.set_ylim((cxlim[0] , cxlim[1] - 1))
        caxis_DI(ax3, 'blue')
        ##PCA Tracks
        ax4 = fig.add_axes([Left + width + 0.1, HB , 0.1 , HH])
        ax4.fill_betweenx(np.arange(PCA.size) , PCA , facecolor = 'green' , edgecolor = 'none' )
        xtick = [PCA.min() + 0.01 , PCA.max() - 0.02]
        ax4.set_xticks(xtick)
        ax4.set_xlabel('PCA',fontsize=20,)
        ax4.set_xlim((PCA.min() , PCA.max()))
        ax4.set_ylim((cxlim[0] , cxlim[1] - 1))
        caxis_DI(ax4, 'green')
        
        #=====================ATAC + RNA Signal================================
        #RNA
        location1 = [Left , HB + HH , width , 0.1]
        color = 'fuchsia'
        Sig_Plot(data = RNA_Data[c],
                 start = startHiC,
                 end = endHiC,
                 chro = g,
                 fig = fig,
                 location = location1,
                 color = color,
                 label = 'RNA')
        
        #ATAC
        location2 = [Left, HB+HH+0.105, width, 0.1]
        color = 'blue'
        Sig_Plot(data = ATAC_Data[c],
                 start = startHiC,
                 end = endHiC,
                 chro = g,
                 fig = fig,
                 location = location2,
                 color = color,
                 label = 'ATAC')
        tmp[c] = fig
    
    
    ccs_fig = tmp['CCS']
    fesc_fig = tmp['fESC']
    nt5_fig = tmp['NT5']
    nt6_fig = tmp['NT6']
    
    figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig]
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
    
        