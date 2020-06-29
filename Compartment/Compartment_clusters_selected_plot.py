# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:53:38 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import cPickle


#---------------------------------------------------Functions-------------------------------------------------------------
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
    ytick = [-0.05 , 0.00 , 0.06]
    ax.set_yticks(ytick)
    ax.set_ylim((-0.06 , 0.07))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)        
    print len(PCA1) , PCA_index1.max()
    return ax
 
def get_union_gene_all(Fil):
    union_type = ({'names':['gene_id','gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3','fESC_1','fESC_2','fESC_3'],
                 'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if ((i['CCS_1'] + i['CCS_2'] +i['CCS_3']) / 3 > 1) or ((i['NT5_1'] + i['NT5_2'] +i['NT5_3'] + i['NT5_4']) / 4 > 1) or ((i['NT6_1'] + i['NT6_2'] +i['NT6_3']) / 3 > 1) or ((i['fESC_1'] + i['fESC_2'] +i['fESC_3']) / 3 > 1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union
    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS','NT5','NT6','fESC']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
R = 200000
res = '200K'



sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

PCFolder = 'H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new'
Outfolder = 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S9_figs'

#interval = [('3' , 32000000 , 40000000)]
genes = {'c2':['Fgf7' , 'Ptgs2' , 'Cxcl5'] , 'c3':['Zeb2' , 'Plxdc2' , 'Sorbs2' ] , 'c4':['Fbn1' , 'Platr3' , 'Mdfic' ] ,
         'c6':['Zscan4c' , 'Fam213a' , 'Rybp-ps'] , 'c7':['Dppa2' , 'Tbx15' , 'Platr15' , 'Tgfbr2'] , 'c8':['Mipol1' , 'Zfp948' , 'Msantd4' , 'Zfp51']}

union_gene_all = get_union_gene_all('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt')
         
interval = {}
for c in genes:
    interval[c] = []
    for gene_name in genes[c]:
        overlap = union_gene_all[union_gene_all['gene_name'] == gene_name][0]
        chro = overlap['chr'].lstrip('chr')
        start = overlap['start']
        end = overlap['end']
        interval[c].append((chro , start - 3000000 , end + 3000000 , gene_name))
        

CCS = np.loadtxt(os.path.join(PCFolder , 'CCS_compartment_200K.txt') , dtype = sig_type)
NT5 = np.loadtxt(os.path.join(PCFolder , 'NT5_compartment_200K.txt') , dtype = sig_type)
NT6 = np.loadtxt(os.path.join(PCFolder , 'NT6_compartment_200K.txt') , dtype = sig_type)
fESC = np.loadtxt(os.path.join(PCFolder , 'fESC_compartment_200K.txt') , dtype = sig_type)



#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 10)
Left1 = 0.2 ; HB = 0.3 ; width1 = 0.45 ; HH = 0.3
Left2 = 0.66 ; width2 = 0.3 

OutFil = 'Abnormal_Compartment_gene_expression.pdf'
pp = PdfPages(os.path.join(Outfolder , OutFil))


for c in ['c2' , 'c3' , 'c4' , 'c6' , 'c7' , 'c8']:
    for i in interval[c]:
        g = i[0]
        start = i[1] // R
        end = i[2] // R
        name = c + '_' + i[3]
        gene = union_gene_all[union_gene_all['gene_name'] == i[3]][0]
        tmp_CCS = CCS[CCS['chr'] == g]
        tmp_NT5 = NT5[NT5['chr'] == g]
        tmp_NT6 = NT6[NT6['chr'] == g]
        tmp_fESC = fESC[fESC['chr'] == g]
        sig1 = tmp_CCS['score'][start:end]
        sig2 = tmp_NT5['score'][start:end]
        sig3 = tmp_NT6['score'][start:end]
        sig4 = tmp_fESC['score'][start:end]
        ##PCA Tracks
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left1, HB , width1 , 0])
        xtick = [i[1] , i[1] + 3000000 , i[2] - 3000000 , i[2]]
        labels = [properU(xtick[0]) , 'S' , 'E' , properU(xtick[-1])]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)

        ax.set_xlabel(name , fontsize=20)
        ax.set_xlim((i[1] , i[2]))
        ax1 = add_ax(sig4 , [Left1, HB , width1 , 0.075] , 'fESC')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig3 , [Left1, HB + 0.075 , width1 , 0.075] , 'NT6')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig2 , [Left1, HB + 0.15, width1 , 0.075] , 'NT5')
        caxis_S_horizontal(ax3)
        ax4 = add_ax(sig1 , [Left1, HB + 0.225, width1 , 0.075] , 'CCS')
        caxis_S_horizontal(ax4)
        
        w = np.log2(np.array([gene['CCS_1'] , gene['CCS_2'] , gene['CCS_3']]) + 1)
        x = np.log2(np.array([gene['NT5_1'] , gene['NT5_2'] , gene['NT5_3'], gene['NT5_4']]) + 1)
        y = np.log2(np.array([gene['NT6_1'] , gene['NT6_2'] , gene['NT6_3']]) + 1)
        z = np.log2(np.array([gene['fESC_1'] , gene['fESC_2'] , gene['fESC_3']]) + 1)
        ax5 = fig.add_axes([Left2, HB , width2 , HH])
        ax5.scatter([1,1,1] , w , s = 700)
        ax5.scatter([2,2,2,2] , x , s = 700)
        ax5.scatter([3,3,3] , y , s = 700)
        ax5.scatter([4,4,4] , z , s = 700)
        ax5.scatter([1,2,3,4] , [w.mean() , x.mean() , y.mean() , z.mean()] , marker = '_' , s = 1000 , c = 'black')
        ax5.set_xlim((0.5,5))
        ax5.set_ylim((0,10))
        ax5.set_xticks([1,2,3,4])
        ax5.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'fESC'] ,fontsize = 20)
#        ax5.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
        
            
        pp.savefig(fig)
        plt.close(fig)
pp.close()
        


