# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 22:52:27 2019

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
import math


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
        if n <17:
            DIs.append(0)
        elif n >= 17 and n <= int(shape[0])-18:
            if len(j[j!=0])/int(shape[0]) < 0.05:
                bias = 0
            else:
                up = j[n-17:n] 
                down = j[n+1:n+18]
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

def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 20)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')
    
def caxis_S(ax):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'x', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 12)
    ax.tick_params(axis = 'y', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 12)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color('black')
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
    s = round(pos / 1000000 , 2)
    
    return str(s) + 'M'

    
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
    ax.set_xlim(0 ,len(sig_data) - 1)
    ax.set_ylim(min(sig_data) , max(sig_data))
    ax.set_ylabel(label,fontsize = 30 ,rotation = 'horizontal',labelpad = 50)
    caxis_S(ax)


def standard_axes_lim(N,*figs):
    """
    """    
    lim_min = []
    lim_max = []
    for s_f in figs:
        lim_min.append(s_f.axes[N].get_ylim()[0])
        lim_max.append(s_f.axes[N].get_ylim()[1])
    
    for s_f in figs:
        s_f.axes[N].set_ylim(min(lim_min),max(lim_max))
        s_f.axes[N].set_yticks([round(min(lim_min)/2 , 1) , 0 , round(max(lim_max)/2 , 1)])
        
        
def diff_sigcolors(sig):
    sig1 = np.zeros(len(sig))
    sig2 = np.zeros(len(sig))
    for x in range(len(sig)):
        if sig[x] > 0:
            sig1[x] = sig[x]
        else:
            sig2[x] = sig[x]
    return sig1 , sig2

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
    loop_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})
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


TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'

chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' ,'X']
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
R = 40000
res = '40K'
#interval = [('17' ,  35506018 , 35510772 , 'Oct4') , ('3' , 34650405 , 34652461 , 'Sox2') , ('4' , 55527143 , 55532466 , 'Klf4') , 
#            ('15' , 61985391 , 61990374 , 'cMyc') , ('6' , 122707489 , 122714633 , 'Nanog') , ('12' , 86361117 , 86521628 , 'Esrrb') , 
#            ('16' , 92601466 , 92826149 , 'Runx1') , ('19' , 5447703 , 5455945 , 'Fra1') , ('7' , 35119293 , 35121928 , 'Cebpa') , 
#            ('2' , 167688915 , 167690418 , 'Cebpb') , ('12' ,  85473890 , 85477273 , 'Fos') , ('4' , 95049034 , 95052222 , 'Jun') , 
#            ('9' , 32636221 , 32757820 , 'Ets1') , ('16' , 95702075 , 95721051 , 'Ets2') , ('7' , 112679318 , 112906807 , 'Tead1') , 
#            ('7' , 45215753 , 45233644 , 'Tead2') , ('17' , 28331671 , 28350805 , 'Tead3') , ('6' , 128224288 , 128300823 , 'Tead4')]

#interval = [('17' ,  35280000 , 36360000 , 'Oct4') , ('3' , 33760000 , 35680000 , 'Sox2') , ('4' , 54680000, 56840000 , 'Klf4') , 
#            ('15' , 59960000 , 63960000 , 'cMyc') , ('6' , 121400000 , 123320000 , 'Nanog') , ('12' , 86080000 , 87120000 , 'Esrrb') , 
#            ('16' , 92320000 , 93720000 , 'Runx1') , ('19' , 5040000 , 6440000 , 'Fra1') , ('5' , 31160000 , 33160000 , 'Fra2')]
#interval = [('1' ,  65000000 , 70000000 ) , ('1' , 180000000 , 186000000), ('2' , 50000000 , 58500000), ('2' , 92500000 , 100000000), ('2' , 103000000 , 110000000), ('2' , 170000000 , 175000000)]
interval = [('2' , 105500000 , 110000000)]
interval_20K = ('2' , 106900000 , 108000000)

size = (12, 12)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
 
union_type = np.dtype({'names':['chr' , 'start' ,'end' ,	'CCS' ,	'NT5'	, 'NT6' , 'fESC'],
                      'formats':['S4' , np.int , np.int , np.float , np.float , np.float , np.float]})
    

    




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

#HiC Data Process
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
#compartment Process
CCS_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/CCS_compartment_200K.txt')
NT5_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/NT5_compartment_200K.txt')
NT6_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/NT6_compartment_200K.txt')
fESC_PCA = pca_To_40K('D:/Workspace_New/data/HiC/Compartment/compartment_new/fESC_compartment_200K.txt')
PCA_Data = {'CCS':CCS_PCA,
            'fESC':fESC_PCA,
            'NT5':NT5_PCA,
            'NT6':NT6_PCA}

#Loop Process

CCS_loop = Get_loops('D:/Workspace_New/data/HiC/Loop/Clustered_CCS_loops_20K.txt')
NT5_loop = Get_loops('D:/Workspace_New/data/HiC/Loop/Clustered_NT5_loops_20K.txt')
NT6_loop = Get_loops('D:/Workspace_New/data/HiC/Loop/Clustered_NT6_loops_20K.txt')
fESC_loop = Get_loops('D:/Workspace_New/data/HiC/Loop/Clustered_fESC_loops_20K.txt')
Loop_Data = {'CCS':CCS_loop,
            'fESC':fESC_loop,
            'NT5':NT5_loop,
            'NT6':NT6_loop} 

union_loop = np.loadtxt('D:/Workspace_New/data/HiC/Loop/loop_strength_point.txt' , dtype = union_type , skiprows = 1)



tmp = {}
pp1 = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig2\\heatmap_40K\\Selected_CCS_40K_Heatmap_stable_window_wise_600K_4.pdf')
pp2 = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig2\\heatmap_40K\\Selected_NT5_40K_Heatmap_stable_window_wise_600K_4.pdf')
pp3 = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig2\\heatmap_40K\\Selected_NT6_40K_Heatmap_stable_window_wise_600K_4.pdf')
pp4 = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig2\\heatmap_40K\\Selected_fESC_40K_Heatmap_stable_window_wise_600K_4.pdf')

for i in interval:
    print i
    g = i[0]
    startHiC = i[1] // 40000
    endHiC = i[2] // 40000
    for c in cell:
        TadFil = c + '_40K_allreps_union_small_domain.txt'
        TadSource = os.path.join(TadFolder , TadFil)
        TadData = np.loadtxt(TadSource , usecols = (0 , 1 , 2) ,dtype = tad_type)
        lib = HiC_Data[c][g]
        DIData = DI_Data[c][g]
        tads = TadData[TadData['chr'] == g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        DI = np.array(DIData[startHiC : endHiC])
        
        
        
        #=============HeatMap + colorbar + DI=================================
        
        fig = plt.figure(figsize = size)
        ax1 = fig.add_axes([Left  , HB , width , HH])
        sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax1.get_xlim()
        cylim = ax1.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 4).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labelsy = [properU(p) for p in pos]
        labelsx = [properU(p) for p in pos[:3]]
        
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labelsx)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labelsy, rotation = 'horizontal' , fontsize = 30)
#        ax1.set_xlabel(c , fontsize = 30 ,labelpad = 30)
        
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
        col = '#000000'
        for corner in pairs:
            if (corner[0] <= 0):
                ax1.plot([0, corner[1]], [corner[1], corner[1]], color = col , linewidth = 5)
                ax1.plot([corner[1], corner[1]], [0, corner[1]], color = col , linewidth = 5)
            elif (corner[1] >= interval):
                ax1.plot([corner[0], corner[0]], [corner[0], interval] , color = col , linewidth = 5)
                ax1.plot([corner[0], interval], [corner[0], corner[0]] , color = col , linewidth = 5)
            else:
                ax1.plot([corner[0], corner[0]], [corner[0], corner[1]] , color = col , linewidth = 5)
                ax1.plot([corner[0], corner[1]], [corner[0], corner[0]] , color = col , linewidth = 5)
                ax1.plot([corner[0], corner[1]], [corner[1], corner[1]] , color = col , linewidth = 5)
                ax1.plot([corner[1], corner[1]], [corner[0], corner[1]] , color = col , linewidth = 5)
        ## Selected interval
        interval_s = interval_20K[1]//R - startHiC 
        interval_e = interval_20K[2]//R - startHiC 
        if g == interval_20K[0]:
            ax1.plot([interval_s, interval_s], [interval_s, interval_e] , color = '#7093DB' , linewidth = 10 , linestyle = '--')
            ax1.plot([interval_s, interval_e], [interval_s, interval_s] , color = '#7093DB' , linewidth = 10 , linestyle = '--')
            ax1.plot([interval_s, interval_e], [interval_e, interval_e] , color = '#7093DB' , linewidth = 10 , linestyle = '--')
            ax1.plot([interval_e, interval_e], [interval_s, interval_e] , color = '#7093DB' , linewidth = 10 , linestyle = '--')
        else:
            pass
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)                    
        caxis_H(ax1)
        ## Colorbar
        ax2 = fig.add_axes([Left + 0.6 , HB - 0.11 , 0.1 , 0.035])
        fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0 , int(vmax)])
        caxis_colorbar(ax2)

        #DI Tracks
        DI_index , DI = UpdateDI(DI)
        ax3 = fig.add_axes([Left , HB + HH , width , 0.075])
        ax3.fill_between(DI_index , DI , where = DI >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax3.fill_between(DI_index , DI , where = DI <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax3.set_xlim(0 , DI_index.max())
#        ax3.set_ylabel('DI',fontsize=30,rotation = 'horizontal' , labelpad = 50)
        
                
        caxis_S_horizontal(ax3,'black')
        tmp[c] = fig
    ccs_fig = tmp['CCS']
    fesc_fig = tmp['fESC']
    nt5_fig = tmp['NT5']
    nt6_fig = tmp['NT6']
    
    figs = [ccs_fig,fesc_fig,nt5_fig,nt6_fig]
    standard_axes_lim(2,*figs)
    
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
