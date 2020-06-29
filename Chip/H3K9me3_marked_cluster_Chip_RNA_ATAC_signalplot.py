# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 17:37:55 2020

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
import seaborn as sns
from scipy.interpolate import  interp1d
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
            end = line['end'] // 100
            for i in range(start,end):
                New_Data[g][i] += line['value']
    
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
    ax = fig.add_axes([Left  , HB + 0.07 * n , width , 0.07])
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
        
        

signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
                    
chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']


Chips = pyBigWig.open('/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw')
Inputs = pyBigWig.open('/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw')
CCS_RNA = pyBigWig.open('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/CCS_RNA_10bp.bw')
NT5_RNA = pyBigWig.open('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/NT5_RNA_10bp.bw')
NT6_RNA = pyBigWig.open('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/NT6_RNA_10bp.bw')
fESC_RNA = pyBigWig.open('/public/home/xxli/data/BDF1_New/RNA/workspace/signal/normalization/bigwig_10bp/fESC_RNA_10bp.bw')
CCS_ATAC = pyBigWig.open('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/CCS_ATAC.bw')
NT5_ATAC = pyBigWig.open('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/NT5_ATAC.bw')
NT6_ATAC = pyBigWig.open('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/NT6_ATAC.bw')
fESC_ATAC = pyBigWig.open('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bigwig_10bp/fESC_ATAC.bw')


pc1 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster1.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc2 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster2.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster3.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc4 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster4.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc5 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster5.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc6 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster6.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc7 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster7.txt' , usecols = (0 , 1, 2) , dtype = pc_type )
pc8 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster8.txt' , usecols = (0 , 1, 2) , dtype = pc_type )




chips = Sig_To_100bp(Chips , 'chr')
inputs = Sig_To_100bp(Inputs , 'chr')
ccs_rna = Sig_To_100bp(CCS_RNA , '')
nt5_rna = Sig_To_100bp(NT5_RNA , '')
nt6_rna = Sig_To_100bp(NT6_RNA , '')
fesc_rna = Sig_To_100bp(fESC_RNA , '')
ccs_atac = Sig_To_100bp(CCS_ATAC , '')
nt5_atac = Sig_To_100bp(NT5_ATAC , '')
nt6_atac = Sig_To_100bp(NT6_ATAC , '')
fesc_atac = Sig_To_100bp(fESC_ATAC , '')


pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }


size = (12, 12)
Left = 0.2 ; HB = 0.15 ; width = 0.7 ; HH = 0.7




for c in ['c2' , 'c3' , 'c4' , 'c6' , 'c7' , 'c8']:
    pp = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/Compartment_' + c + '_H3K9me3_RNA_ATAC_signals.pdf')
    for i in pc_data[c]:
        chro = i['chr']
        start = i['start'] // 100
        end = i['end'] // 100
        sig1 = chips[chro][start:end]
        sig2 = inputs[chro][start:end]
        sig3 = ccs_rna[chro][start:end]
        sig4 = nt5_rna[chro][start:end]
        sig5 = nt6_rna[chro][start:end]
        sig6 = fesc_rna[chro][start:end]
        sig7 = ccs_atac[chro][start:end]
        sig8 = nt5_atac[chro][start:end]
        sig9 = nt6_atac[chro][start:end]
        sig10 = fesc_atac[chro][start:end]
        rna_max = max([sig3.max() , sig4.max() , sig5.max() , sig6.max()]) + 5
        atac_max = max([sig7.max() , sig8.max() , sig9.max() , sig10.max()]) + 5
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , 0])
        xtick = [0 , len(sig1)]
        labels = [properU(start * 100) , properU(end * 100)]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)
        ax.set_xlabel('chr' + chro)
        ax1 = Get_signal_ax(fig , 0 , sig1 , 'green' , 'K9_chip' , sig1.max())
        ax2 = Get_signal_ax(fig , 1 , sig2 , 'green' , 'K9_input' , sig1.max())
        ax3 = Get_signal_ax(fig , 2 , sig3 , 'crimson' , 'CCS_RNA' , rna_max)
        ax4 = Get_signal_ax(fig , 3 , sig4 , 'crimson' , 'NT5_RNA' , rna_max)
        ax5 = Get_signal_ax(fig , 4 , sig5 , 'crimson' , 'NT6_RNA' , rna_max)
        ax6 = Get_signal_ax(fig , 5 , sig6 , 'crimson' , 'fESC_RNA' , rna_max)
        ax7 = Get_signal_ax(fig , 6 , sig7 , 'blue' , 'CCS_ATAC' , atac_max)
        ax8 = Get_signal_ax(fig , 7 , sig8 , 'blue' , 'NT5_ATAC' , atac_max)
        ax9 = Get_signal_ax(fig , 8 , sig9 , 'blue' , 'NT6_ATAC' , atac_max)
        ax10 = Get_signal_ax(fig , 9 , sig10 , 'blue' , 'fESC_ATAC' , atac_max)
        pp.savefig(fig)
    pp.close()
        
        

        
        
        
