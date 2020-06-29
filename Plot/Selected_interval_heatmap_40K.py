# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:23:05 2018

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
                
                
#-----------------------------------------------------Functions---------------------------------------------------------
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
    

#--------------------------------------------------Files----------------------------------------------------------------


TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\stable_window_wsize_800K'
ATACFolder = 'D:\\Workspace_New\\data\\ATAC\\signal\\normalization\\bedgraph_1K'
RNAFolder = 'D:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K'
OutFolder = 'D:\\Workspace_New\\Plot\\Figure2\\somatic_gene'

chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' ,'X']
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
R = 40000
res = '40K'
#interval = [('17' ,  35506018 , 35510772 , 'Oct4') , ('3' , 34650405 , 34652461 , 'Sox2') , ('4' , 55527143 , 55532466 , 'Klf4') , 
#            ('15' , 61985391 , 61990374 , 'cMyc') , ('6' , 122707489 , 122714633 , 'Nanog') , ('12' , 86361117 , 86521628 , 'Esrrb') , 
#            ('16' , 92601466 , 92826149 , 'Runx1') , ('19' , 5447703 , 5455945 , 'Fra1') , ('7' , 35119293 , 35121928 , 'Cebpa') , 
#            ('2' , 167688915 , 167690418 , 'Cebpb')]

interval = [('12' ,  85473890 , 85477273 , 'Fos') , ('4' , 95049034 , 95052222 , 'Jun') , ('9' , 32636221 , 32757820 , 'Ets1') , 
            ('16' , 95702075 , 95721051 , 'Ets2') , ('7' , 112679318 , 112906807 , 'Ted1') , ('7' , 45215753 , 45233644 , 'Ted2') , 
            ('17' , 28331671 , 28350805 , 'Tead3') , ('6' , 128224288 , 128300823 , 'Tead4')]
size = (12, 12)
Left = 0.15 ; HB = 0.15 ; width = 0.6 ; HH = 0.6

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
    
    
 
for c in cell:
    OutFil = 'Selected_' + c + '_40K_Heatmap_stable_window_wise_800K.pdf'
    pp = PdfPages(os.path.join(OutFolder , OutFil))
    HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\' + c + '_40K'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource)
    TadFil = c + '_40K_allreps.txt'
#    DIFil = c + '_40K_DI.txt'
    TadSource = os.path.join(TadFolder , TadFil)
#    DISource = os.path.join(TadFolder , DIFil)
    TadData = np.loadtxt(TadSource , usecols = (0 , 1 , 2) ,dtype = tad_type)
    ATACFil = c + '_ATAC_1K.bedgraph'
    RNAFil = c + '_RNA_1K.bedgraph'
    ATACSource = os.path.join(ATACFolder , ATACFil)
    RNASource = os.path.join(RNAFolder , RNAFil)
    ATACData = np.loadtxt(ATACSource , dtype=sig_type)
    RNAData = np.loadtxt(RNASource , dtype=sig_type)
    for i in interval:
        g = i[0]
        startHiC = i[1] // R - 25 
        endHiC = i[2] // R + 25 
        lib = HiCData['Matrix'][()][g]
        DIData = CalDI(HiCData['Matrix'][()][g])
        tads = TadData[TadData['chr'] == g]
        matrix = lib[startHiC:endHiC , startHiC:endHiC]
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        DI = np.array(DIData[startHiC : endHiC])
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [((startHiC + t) * R) for t in ticks]
        labels = [properU(p) for p in pos[:4]]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, rotation = 'horizontal')
        ax.set_xlabel(i[-1] , fontsize = 30 )
        
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
                ax.plot([0, corner[1]], [corner[1], corner[1]], color = '#000000' , linewidth = 1)
                ax.plot([corner[1], corner[1]], [0, corner[1]], color = '#000000' , linewidth = 1)
            elif (corner[1] >= interval):
                ax.plot([corner[0], corner[0]], [corner[0], interval] , color = '#000000' , linewidth = 1)
                ax.plot([corner[0], interval], [corner[0], corner[0]] , color = '#000000' , linewidth = 1)
            else:
                ax.plot([corner[0], corner[0]], [corner[0], corner[1]] , color = '#000000' , linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[0], corner[0]] , color = '#000000' , linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[1], corner[1]] , color = '#000000' , linewidth = 1)
                ax.plot([corner[1], corner[1]], [corner[0], corner[1]] , color = '#000000' , linewidth = 1)
                        
#b3b3b3
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        caxis_H(ax)
        ## Colorbar
        ax = fig.add_axes([Left + 0.5 , HB - 0.8 , 0.1 , 0.035])
        cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
        cbar.set_ticks([0 , int(vmax/2) , int(vmax)])

        ##DI Tracks
        ax = fig.add_axes([Left + width , HB , 0.1 , HH])
        ax.fill_betweenx(np.arange(DI.size) , DI , facecolor = 'blue' , edgecolor = 'none' )
        xtick = [DI.min() + 1 , DI.max() - 1]
        ax.set_xticks(xtick)
        ax.set_xlabel('DI',fontsize=20)
        ax.set_xlim((DI.min() , DI.max()))
        ax.set_ylim((cxlim[0] , cxlim[1] - 1))
        caxis_DI(ax, 'blue')
        ##Signal Tracks
        atac = ATACData[ATACData['chr'] == g]
#        mask2 = ((i[1] // 1000) <= atac['start'] // 1000) & ((i[2] // 1000) >= atac['end'] // 1000)
#        atac = atac[mask2]
        rna = RNAData[RNAData['chr'] == g]
#        mask3 = ((i[1] // 1000) <= rna['start'] // 1000) & ((i[2] // 1000) >= rna['end'] // 1000)
#        rna = rna[mask3]
        signal = {'ATAC':atac , 'RNA':rna}
        step = 0
        for s in ['RNA','ATAC']:
            length = mm[g]//1000
            new = np.zeros(length)
            for ss in signal[s]:
                start = ss['start'] // 1000 
                end = ss['end'] // 1000
                score = ss['score']
                n = end - start
                for nn in np.arange(n):
                    new[start + nn] = score
            new = new[(startHiC * R) // 1000 : (endHiC * R) // 1000]
            if i[-1] == 'Fos':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 200
                else:
                    color = 'blue'
                    vm = 25
            elif i[-1] == 'Jun':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 100
                else:
                    color = 'blue'
                    vm = 30
            elif i[-1] == 'Ets1':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 30
                else:
                    color = 'blue'
                    vm = 15
            elif i[-1] == 'Ets2':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 200
                else:
                    color = 'blue'
                    vm = 25
            elif i[-1] == 'Tead1':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 75
                else:
                    color = 'blue'
                    vm = 25
            elif i[-1] == 'Tead2':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 200
                else:
                    color = 'blue'
                    vm = 30
            elif i[-1] == 'Tead3':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 500
                else:
                    color = 'blue'
                    vm = 20
            elif i[-1] == 'Tead4':
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = 200
                else:
                    color = 'blue'
                    vm = 20
    
            else:
                if s == 'RNA':
                    color = 'fuchsia'
                    vm = new.max()
                else:
                    color = 'blue'
                    vm = new.max()
                    
            
            ax = fig.add_axes([Left , HB + HH + step , width , 0.1])
            ax.fill_between(np.arange(len(new)) , new , facecolor = color , edgecolor = 'none' )
            yticks = [0 , int(vm * 0.45) , int(vm * 0.85)]
            ax.set_yticks(yticks)
            ax.set_ylabel(s,fontsize=15)
            ax.set_ylim((0 , vm))
            ax.set_xlim((0 , len(new) - 1))
            caxis_S(ax, color)
            step += 0.1
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
            
#fig = plt.figure(figsize = size)
#ax = fig.add_axes([Left  , HB , width , HH])
#ax.imshow(matrix)
#ax.plot([0,0] , [1,2])       
#caxis_H(ax)