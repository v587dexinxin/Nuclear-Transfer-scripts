# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:02:30 2020

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

#---------------------------------------------------Functions-------------------------------------------------------------
def Get_PC(PC_fil):
    '''
    '''
    data_type = np.dtype({'names':['chr' , 'start' , 'end'],
                          'formats':['U8' , np.int , np.int]})
    data = np.loadtxt(PC_fil , skiprows=1 , usecols=(0 , 1 , 2) , dtype = data_type)
    return data

    
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
    
    
def caxis_S_vertical(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'y', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
def caxis_S_horizontal(ax):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

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

def add_ax(sig , loc , cell):
    PCA_index1 , PCA1 = UpdateDI(sig)
    ax = fig.add_axes(loc)
    ax.fill_between(PCA_index1 , PCA1 , where = PCA1 >= 0 , facecolor = 'gold' , edgecolor = 'none' )
    ax.fill_between(PCA_index1 , PCA1 , where = PCA1 <= 0 , facecolor = 'midnightblue' , edgecolor = 'none' )
    ax.set_xlim((0 , PCA_index1.max()))
    ytick = [-0.04 , 0.00 , 0.04]
    ax.set_yticks(ytick)
    ax.set_ylim((-0.1 , 0.1))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)        
    print (len(PCA1) , PCA_index1.max())
    return ax

def add_ax(sig , loc , cell):
    sig1 = [] ; sig2 = [] ; index1 = [] ; index2 = []
    for i in range(len(sig)):
        if sig[i] > 0:
            index1.append(i)
            sig1.append(sig[i])
        elif sig[i] < 0:
            index2.append(i)
            sig2.append(sig[i])
        else:
            pass
    index1 = np.array(index1)
    index2 = np.array(index2)
    sig1 = np.array(sig1)
    sig2 = np.array(sig2)

    
    ax = fig.add_axes(loc)
    ax.bar(index1 , sig1 , 1 , color = 'gold')
    ax.bar(index2 , sig2 , 1 , color = 'midnightblue')
    # ax.fill_between(PCA_index1 , PCA1 , where = PCA1 >= 0 , facecolor = 'gold' , edgecolor = 'none' )
    # ax.fill_between(PCA_index1 , PCA1 , where = PCA1 <= 0 , facecolor = 'midnightblue' , edgecolor = 'none' )
    ax.set_xlim((-0.5 , len(sig) - 0.5))
    ytick = [-0.04 , 0.00 , 0.04]
    ax.set_yticks(ytick)
    ax.set_ylim((-0.1 , 0.1))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)        
    print (len(sig))
    return ax

        
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'F35':3 , 'F40':4}
cell = ['CCS','NT5','NT6','F35','F40']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
R = 200000
res = '200K'

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U8',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['U8',np.int , np.int , np.float]})

PCFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new'
Outfolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\plot\\Abnormal_compartment_plot'

#interval = [('3' , 32000000 , 40000000)]
# interval = [('5' , 88000000 , 94000000)]
#interval = [('16' , 5000000 , 47000000) , ('18' , 45000000 , 90000000)]

CCS = np.loadtxt(os.path.join(PCFolder , 'CCS_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT5 = np.loadtxt(os.path.join(PCFolder , 'NT5_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT6 = np.loadtxt(os.path.join(PCFolder , 'NT6_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F35 = np.loadtxt(os.path.join(PCFolder , 'F35_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F40 = np.loadtxt(os.path.join(PCFolder , 'F40_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)

NTs = np.loadtxt(os.path.join(PCFolder , 'NTs_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
fESC = np.loadtxt(os.path.join(PCFolder , 'fESC_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)






#--------------------------------------------compartment_classify-----------------------------------------------


c1 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster1_pc1.txt')
c2 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster2_pc1.txt')
c3 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster3_pc1.txt')
c4 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster4_pc1.txt')
c5 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster5_pc1.txt')
c6 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster6_pc1.txt')
c7 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster7_pc1.txt')
c8 = Get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster8_pc1.txt')


cl = {'A_To_B_Reprogramming_compartment' : c1 , 'A_To_B_Partial_compartment' : c2 , 
      'A_To_B_Over_compartment' : c3 , 'A_To_B_Resistant_compartment' : c4 ,
      'B_To_A_Reprogramming_compartment' : c5 , 'B_To_A_Partial_compartment' : c6 ,
      'B_To_A_Over_compartment' : c7 , 'B_To_A_Resistant_compartment' : c8 }


#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 10)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7

for c in list(cl.keys())[:2]:
    OutFil = c + '_200K_1.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    
    
    for i in cl[c]:
        g = i['chr']
        start = i['start'] // R - 20
        if start < 0:
            start = 0
        end = i['end'] // R + 20
        tmp_CCS = CCS[CCS['chr'] == g]
        tmp_NT5 = NT5[NT5['chr'] == g]
        tmp_NT6 = NT6[NT6['chr'] == g]
        tmp_F35 = F35[F35['chr'] == g]
        tmp_F40 = F40[F40['chr'] == g]
        sig1 = tmp_CCS['score'][start:end]
        sig2 = tmp_NT5['score'][start:end]
        sig3 = tmp_NT6['score'][start:end]
        sig4 = tmp_F35['score'][start:end]
        sig5 = tmp_F40['score'][start:end]
        ##PCA Tracks
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left, HB , width , 0])
        
        for spine in ['right', 'top']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                       right = False, labelbottom = True, labeltop = False,
                       labelleft = False, labelright = False , labelsize = 23)
        xtick = [0 , i['start'] // R - start + 1 , len(sig4)]
        pos = [((start + t) * R) for t in xtick]
        labels = [properU(p) for p in pos]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)
        ax.set_xlabel(g , fontsize=20)
        ax.set_xlim((0 , len(sig5)))
        ax1 = add_ax(sig5 , [Left, HB , width , 0.1] , 'F40')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig4 , [Left, HB + 0.1 , width , 0.1] , 'F35')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig3 , [Left, HB + 0.2, width , 0.1] , 'NT6')
        caxis_S_horizontal(ax3)
        ax4 = add_ax(sig2 , [Left, HB + 0.3, width , 0.1] , 'NT5')
        caxis_S_horizontal(ax4)
        ax5 = add_ax(sig1 , [Left, HB + 0.4, width , 0.1] , 'CCS')
        caxis_S_horizontal(ax5)
            
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
        
    
#-----------------------------------------Plot_NTs_fESC_merged---------------------------------------------------------------------------

size = (12, 10)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7

for c in ['B_To_A_Partial_compartment']:
    OutFil = c + '_200K_NTs_fESC_merged.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    
    
    for i in cl[c]:
        g = i['chr']
        start = i['start'] // R - 20
        if start < 0:
            start = 0
        end = i['end'] // R + 20
        tmp_CCS = CCS[CCS['chr'] == g]
        tmp_NTs = NTs[NTs['chr'] == g]
        tmp_fESC = fESC[fESC['chr'] == g]
        sig1 = tmp_CCS['score'][start:end]
        sig2 = tmp_NTs['score'][start:end]
        sig3 = tmp_fESC['score'][start:end]

        ##PCA Tracks
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left, HB , width , 0])
        
        for spine in ['right', 'top']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                       right = False, labelbottom = True, labeltop = False,
                       labelleft = False, labelright = False , labelsize = 23)
        xtick = [0 , i['start'] // R - start + 1 , len(sig1)]
        pos = [((start + t) * R) for t in xtick]
        labels = [properU(p) for p in pos]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)
        ax.set_xlabel(g , fontsize=20)
        ax.set_xlim((0 , len(sig1)))
        ax1 = add_ax(sig3 , [Left, HB , width , 0.1] , 'fESC')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig2 , [Left, HB + 0.1, width , 0.1] , 'NTs')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig1 , [Left, HB + 0.2, width , 0.1] , 'CCS')
        caxis_S_horizontal(ax3)
            
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
           
