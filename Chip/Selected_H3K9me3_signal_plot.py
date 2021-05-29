# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:36:40 2021

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
import pyBigWig
import seaborn as sns
from scipy.interpolate import  interp1d
from scipy import stats

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


def Sig_To_100bp(signal , chro):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals(chro + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            # end = line['end'] // 100
            # for i in range(start,end):
            New_Data[g][start] += line['value']
    
    return New_Data


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


    
def Get_signal_ax(fig , n , sig , color , labels , ylim):
    ax = fig.add_axes([Left  , HB + 0.05 * n , width , 0.05])
    ax.fill_between(np.arange(len(sig)),sig, facecolor = color, edgecolor = 'none')
    ax.set_xlim((0 , len(sig)))
    ax.set_ylim((0 , ylim))
    ax.set_ylabel(labels,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)
    return ax

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
        
def pca_To_200K(fil):
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
            New_Data[c].extend([i['PCA']])
            
    
    return New_Data


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


def Get_PCA_ax(fig , n , PCA , color , labels):
    ##PCA Tracks
    # PCA_index , PCA = UpdateDI(PCA)
    x = [1 + w * 0.5 for w in range(len(PCA))]
    y = PCA
    y1 = np.zeros(len(PCA))
    y2 = np.zeros(len(PCA))
    for i in range(len(y)):
        if y[i] > 0:
            y1[i] = y[i]
        else:
            y2[i] = y[i]
    ax = fig.add_axes([Left  , HB + 0.07 * n , width , 0.07])
    ax.bar(x , y1 , 0.5 , color = '#E47833' , edgecolor = 'none' )
    ax.bar(x , y2 , 0.5 , color = '#7093DB' , edgecolor = 'none' )
    
    ax.set_xlim(0.5 , max(x)+0.5)
    ax.set_ylim(-0.06 , 0.02)
    ax.set_yticks([-0.03 , 0 , 0.015])
    ax.set_ylabel(labels + '_PC1')
    caxis_S_horizontal(ax,color)


signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
tad_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                       'formats':['S8' , np.int , np.int]})
                    
chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']


CCS_PCA = pca_To_200K('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/CCS_Traditonal_PC_200K_Compartment_200K.txt')
NT5_PCA = pca_To_200K('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/NT5_Traditonal_PC_200K_Compartment_200K.txt')
NT6_PCA = pca_To_200K('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/NT6_Traditonal_PC_200K_Compartment_200K.txt')
F35_PCA = pca_To_200K('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/F35_Traditonal_PC_200K_Compartment_200K.txt')        
F40_PCA = pca_To_200K('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/F40_Traditonal_PC_200K_Compartment_200K.txt') 




CCS_Chips = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_chip_50bp.bw')
CCS_input = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_input_50bp.bw')

CCS_Chips = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw')
CCS_input = pyBigWig.open('/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw')

ESC_Chips = pyBigWig.open('/public/home/lixinxin/data/literature/cell_cooperative/mapping/Unique_ESC_H3K9me3_50bp.bw')
ESC_input = pyBigWig.open('/public/home/lixinxin/data/literature/cell_cooperative/mapping/Unique_ESC_H3K9me3_Input_50bp.bw')

ccs_chips = Sig_To_100bp(CCS_Chips , 'chr')
ccs_input = Sig_To_100bp(CCS_input , 'chr')
esc_chips = Sig_To_100bp(ESC_Chips , '')
esc_input = Sig_To_100bp(ESC_input , '')




CCS_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/CCS_Speci_TADs.txt' , dtype = tad_type)
fESC_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/HiTAD/HiTAD_new/xtwang/bottom_domain/TAD_classify_weix/fESC_Speci_TADs.txt' , dtype = tad_type)

H3K9me3_CCS_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_TADs/CCS_Speci_CCS_H3K9me3_marked_TADs_fc1.4.txt' , dtype = tad_type)
H3K9me3_fESC_Speci_Tad = np.loadtxt('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_TADs/fESC_Speci_fESC_H3K9me3_marked_TADs_fc1.4.txt' , dtype = tad_type)

CCS_unmarked = []
fESC_unmarked = []

for i in CCS_Speci_Tad:
    tmp = H3K9me3_CCS_Speci_Tad[H3K9me3_CCS_Speci_Tad['chr'] == i[0]]
    mask = (tmp['start'] == i['start'] ) & (tmp['end'] == i['end'] )
    overlap = tmp[mask]
    if overlap.size == 0:
        CCS_unmarked.append(i)
CCS_unmarked = np.array(CCS_unmarked , dtype = CCS_Speci_Tad.dtype)


for i in fESC_Speci_Tad:
    tmp = H3K9me3_fESC_Speci_Tad [H3K9me3_fESC_Speci_Tad ['chr'] == i[0]]
    mask = (tmp['start'] == i['start'] ) & (tmp['end'] == i['end'] )
    overlap = tmp[mask]
    if overlap.size == 0:
        fESC_unmarked.append(i)
fESC_unmarked = np.array(fESC_unmarked , dtype = CCS_Speci_Tad.dtype)



TAD = {'CCS_Speci': H3K9me3_CCS_Speci_Tad , 'fESC_Speci':H3K9me3_fESC_Speci_Tad , 'CCS_Speci_unmarked':CCS_unmarked , 'fESC_Speci_unmaked':fESC_unmarked}




size = (12, 12)
Left = 0.2 ; HB = 0.1 ; width = 0.7 ; HH = 0.7

for c in TAD:
    pp = PdfPages('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_TADs/' + c + '_TADs_H3K9me3_signals_plot_2.pdf')

    for i in TAD[c]:
        chro = i[0]
        start = i[1] // 100
        end = i[2] // 100
        sig1 = ccs_chips[chro][start:end]
        sig2 = ccs_input[chro][start:end]
        sig3 = esc_chips[chro][start:end]
        sig4 = esc_input[chro][start:end]

      
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , 0])
        xtick = [0 , len(sig1)]
        labels = [properU(start * 100) , properU(end * 100)]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)
        ax.set_xlabel('chr' + chro)
        ax1 = Get_signal_ax(fig , 0 , sig4 , 'green' , 'ESC_H3K9me3_Input' , 25)
        ax2 = Get_signal_ax(fig , 1 , sig3 , 'green' , 'ESC_H3K9me3_Chip' , 25)
        ax3 = Get_signal_ax(fig , 2 , sig2 , 'crimson' , 'CCS_H3K9me3_Input' , 25)
        ax4 = Get_signal_ax(fig , 3 , sig1 , 'crimson' , 'CCS_H3K9me3_Chip' , 25)

        
        # ax = fig.add_axes([Left  , HB + 0.07 * 6 , width , 0])
        # xtick = [0 , len(sig1)]
        # labels = [properU(i[1]) , properU(i[2])]
        # ax.set_xticks(xtick)
        # ax.set_xticklabels(labels)
        # ax.set_xlabel('Chr' + chro)
        
        
        pp.savefig(fig)
        plt.close()
    pp.close()
    
    