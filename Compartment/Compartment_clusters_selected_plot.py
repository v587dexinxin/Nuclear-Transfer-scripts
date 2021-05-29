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
# import cPickle


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
    print (len(PCA1) , PCA_index1.max())
    return ax
 

def get_union_gene_all(Fil):
    union_type = ({'names':['gene_id' , 'gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3','F35_1','F35_2','F35_3','F40_1','F40_2','F40_3'],
                 'formats':['U64' , 'U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float, np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if ((i['CCS_1'] + i['CCS_2'] +i['CCS_3']) / 3 > 0.1) or ((i['NT5_1'] + i['NT5_2'] +i['NT5_3']) / 3 > 0.1) or ((i['NT6_1'] + i['NT6_2'] +i['NT6_3']) / 3 > 0.1) or ((i['F35_1'] + i['F35_2'] +i['F35_3']) / 3 > 0.1) or ((i['F40_1'] + i['F40_2'] +i['F40_3']) / 3 > 0.1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union

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
    ytick = [-0.05 , 0.00 , 0.06]
    ax.set_yticks(ytick)
    ax.set_ylim((-0.1 , 0.1))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)        
    print (len(sig))
    return ax

    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'f35':3 , 'f40':4}
cell = ['CCS','NT5','NT6','F35','F40']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
R = 200000
res = '200K'



sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['U4',np.int , np.int , np.float]})

PCFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new'
Outfolder = 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S7'

#interval = [('3' , 32000000 , 40000000)]
genes = {'c1':['Shisa6'] ,'c2':['Cd34' , 'Ica1'] , 'c3':['Zeb2' , 'Gm27761' , 'Zeb2os' , 'Mir5129' , 'Senp6' , 'Gm6081'] , 
         'c4':['Rcan2', 'Acyp2' , 'Dab2', 'Rftn1', 'Fam13a', 'Maf', 'Bmp2k', 'Far2','Btg1' , 'Acyp2' , 'Zc2hc1a' , 'Pop4' , 'Tfg' , 'Nap1l5' , 'Rcan2' , 'Me1' , 'Mmadhc' , 'Gm38230'] , 
         'c5':['Tbx3'] , 
         'c6':['n-R5s45'] , 
         'c7':['Msc' , 'Ildr1' , 'Amph' , 'Schip1' , 'Gm21949' , 'Eri2' , 'Thumpd1' , 'Cyp2b23' , 'Syce1' , 'Zfp951' ,
        'Dppa2' , 'Cnksr3' , 'Mrps11' , 'Mrpl46' , 'Ppa2' , 'AC164424.2' , 'Nup37' , 'Zfp518b' , 'Rexo5' , 'Dcun1d3' , 'Rex2' , 
        'Fam96b' , '6820431F20Rik' , 'Ccdc58' , 'Mme' , 'Nsfl1c' , 'Wdr3' , 'Mier2' , 'Cd38' , 'Fbxo30' , 'Trim41' , 'Msrb2' , 
        'Gm14322' , 'Rrs1'] , 'c8':['Sftpd','Lyst','Cdhr1' , 'AI838599' , 'Shcbp1' , 'Zfp958' , 'Chek1' , 'Ei24','Cdon','Hyls1','Ttk','Gm42959']}


# genes = {'c6':['n-R5s45', 'Cadps', 'AC156025.2', 'Gm38292', 'Gm37474', 'Gm37082',
#        'Gm22060', 'Gm27201', 'Fli1', 'Gm38346', 'Gm43294', 'Hepacam2',
#        'Vps50', 'Gm24065', 'Gpr137b-ps', 'Gm2399', 'Prl2c3', 'CT030166.1',
#        'CT030166.6', 'CT030166.5', 'CT030166.4', 'CT030166.3',
#        'CT030166.2']}

union_gene_all = get_union_gene_all('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_genes.txt')
         
interval = {}
for c in genes:
    interval[c] = []
    for gene_name in genes[c]:
        overlap = union_gene_all[union_gene_all['gene_name'] == gene_name][0]
        chro = overlap['chr'].lstrip('chr')
        start = overlap['start']
        end = overlap['end']
        interval[c].append((chro , start - 3000000 , end + 3000000 , gene_name))
        

CCS = np.loadtxt(os.path.join(PCFolder , 'CCS_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT5 = np.loadtxt(os.path.join(PCFolder , 'NT5_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT6 = np.loadtxt(os.path.join(PCFolder , 'NT6_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F35 = np.loadtxt(os.path.join(PCFolder , 'F35_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F40 = np.loadtxt(os.path.join(PCFolder , 'F40_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)


#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 10)
Left1 = 0.2 ; HB = 0.3 ; width1 = 0.45 ; HH = 0.3
Left2 = 0.66 ; width2 = 0.3 

OutFil = 'Compartment_gene_expression_3.pdf'
pp = PdfPages(os.path.join(Outfolder , OutFil))


for c in ['c1','c2','c3','c4','c5','c6','c7','c8']:
# for c in ['c3']:
    for i in interval[c]:
        g = i[0]
        start = i[1] // R
        end = i[2] // R
        name = c + '_' + i[3] + '_Chr' + g
        gene = union_gene_all[union_gene_all['gene_name'] == i[3]][0]
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
        ax = fig.add_axes([Left1, HB , width1 , 0])
        xtick = [start * R , i[1] + 3000000 , i[2] - 3000000 , end * R]
        labels = [properU(xtick[0]) , 'S' , 'E' , properU(xtick[-1])]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)

        ax.set_xlabel(name , fontsize=20)
        ax.set_xlim((start * R , end * R ))
        ax0 = add_ax(sig5 , [Left1, HB , width1 , 0.075] , 'F40')
        caxis_S_horizontal(ax0)
        ax1 = add_ax(sig4 , [Left1, HB + 0.075 , width1 , 0.075] , 'F35')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig3 , [Left1, HB + 0.15 , width1 , 0.075] , 'NT6')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig2 , [Left1, HB + 0.225, width1 , 0.075] , 'NT5')
        caxis_S_horizontal(ax3)
        ax4 = add_ax(sig1 , [Left1, HB + 0.3, width1 , 0.075] , 'CCS')
        caxis_S_horizontal(ax4)
        
        v = np.log2(np.array([gene['CCS_1'] , gene['CCS_2'] , gene['CCS_3']]) + 1)
        w = np.log2(np.array([gene['NT5_1'] , gene['NT5_2'] , gene['NT5_3']]) + 1)
        x = np.log2(np.array([gene['NT6_1'] , gene['NT6_2'] , gene['NT6_3']]) + 1)
        y = np.log2(np.array([gene['F35_1'] , gene['F35_2'] , gene['F35_3']]) + 1)
        z = np.log2(np.array([gene['F40_1'] , gene['F40_2'] , gene['F40_3']]) + 1)
        ax5 = fig.add_axes([Left2, HB , width2 , HH + 0.075])
        ax5.scatter([1,1,1] , v , s = 700)
        ax5.scatter([2,2,2] , w , s = 700)
        ax5.scatter([3,3,3] , x , s = 700)
        ax5.scatter([4,4,4] , y , s = 700)
        ax5.scatter([5,5,5] , z , s = 700)
        ax5.scatter([1,2,3,4,5] , [v.mean() , w.mean() , x.mean() , y.mean() , z.mean()] , marker = '_' , s = 1000 , c = 'black')
        ax5.set_xlim((0.5,6))
        ax5.set_ylim((0,10))
        ax5.set_xticks([1,2,3,4,5])
        ax5.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40'] ,fontsize = 20)
#        ax5.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
        
            
        pp.savefig(fig)
        plt.close(fig)
pp.close()
        


#----------------------------------------------Plot_Msc---------------------------------------------------------------------------
size = (12, 10)
Left1 = 0.2 ; HB = 0.3 ; width1 = 0.45 ; HH = 0.3
Left2 = 0.66 ; width2 = 0.3 

Outfolder = 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\Fig4'
OutFil = 'Abnormal_Compartment_gene_expression_1.pdf'
pp = PdfPages(os.path.join(Outfolder , OutFil))


for c in ['c7']:
# for c in ['c3']:
    for i in interval[c]:
        g = i[0]
        start = i[1] // R
        end = i[2] // R
        name = c + '_' + i[3] + '_Chr' + g
        gene = union_gene_all[union_gene_all['gene_name'] == i[3]][0]
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
        ax = fig.add_axes([Left1, HB , width1 , 0])
        xtick = [i[1] , 14000000 , 15000000 , 15600000 , 16800000 ,  i[1] + 3000000 , i[2] - 3000000 , i[2]]
        labels = [properU(xtick[0]) , 'o1' , '' , 'o2' , 'o3' , 'S' , 'E' , properU(xtick[-1])]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)

        ax.set_xlabel(name , fontsize=20)
        ax.set_xlim((i[1] , i[2]))
        ax0 = add_ax(sig5 , [Left1, HB , width1 , 0.075] , 'F40')
        caxis_S_horizontal(ax0)
        ax1 = add_ax(sig4 , [Left1, HB + 0.075 , width1 , 0.075] , 'F35')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig3 , [Left1, HB + 0.15 , width1 , 0.075] , 'NT6')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig2 , [Left1, HB + 0.225, width1 , 0.075] , 'NT5')
        caxis_S_horizontal(ax3)
        ax4 = add_ax(sig1 , [Left1, HB + 0.3, width1 , 0.075] , 'CCS')
        caxis_S_horizontal(ax4)
        
        v = np.log2(np.array([gene['CCS_1'] , gene['CCS_2'] , gene['CCS_3']]) + 1)
        w = np.log2(np.array([gene['NT5_1'] , gene['NT5_2'] , gene['NT5_3']]) + 1)
        x = np.log2(np.array([gene['NT6_1'] , gene['NT6_2'] , gene['NT6_3']]) + 1)
        y = np.log2(np.array([gene['F35_1'] , gene['F35_2'] , gene['F35_3']]) + 1)
        z = np.log2(np.array([gene['F40_1'] , gene['F40_2'] , gene['F40_3']]) + 1)
        ax5 = fig.add_axes([Left2, HB , width2 , HH + 0.075])
        ax5.scatter([1,1,1] , v , s = 700)
        ax5.scatter([2,2,2] , w , s = 700)
        ax5.scatter([3,3,3] , x , s = 700)
        ax5.scatter([4,4,4] , y , s = 700)
        ax5.scatter([5,5,5] , z , s = 700)
        ax5.scatter([1,2,3,4,5] , [v.mean() , w.mean() , x.mean() , y.mean() , z.mean()] , marker = '_' , s = 1000 , c = 'black')
        ax5.set_xlim((0.5,6))
        ax5.set_ylim((0,10))
        ax5.set_xticks([1,2,3,4,5])
        ax5.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40'] ,fontsize = 20)
#        ax5.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
        
            
        pp.savefig(fig)
        plt.close(fig)
pp.close()
        


#----------------------------------------------Plot_all---------------------------------------------------------------------------

c2_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\A_To_B_Partial_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)
c3_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\A_To_B_Over_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)
c4_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\A_To_B_Resistant_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)
c6_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\B_To_A_Partial_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)
c7_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\B_To_A_Over_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)
c8_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\B_To_A_Resistant_compartment_related_genes_nofiltering.txt' , skiprows = 1 , dtype = 'U64' , usecols = 1)

genes = {'A_To_B_Partial':c2_gene ,'A_To_B_Over':c3_gene , 'A_To_B_Resistant':c4_gene , 
         'B_To_A_Partial': c6_gene, 'B_To_A_Over': c7_gene, 'B_To_A_Resistant': c8_gene}


# genes = {'c6':['n-R5s45', 'Cadps', 'AC156025.2', 'Gm38292', 'Gm37474', 'Gm37082',
#        'Gm22060', 'Gm27201', 'Fli1', 'Gm38346', 'Gm43294', 'Hepacam2',
#        'Vps50', 'Gm24065', 'Gpr137b-ps', 'Gm2399', 'Prl2c3', 'CT030166.1',
#        'CT030166.6', 'CT030166.5', 'CT030166.4', 'CT030166.3',
#        'CT030166.2']}

union_gene_all = get_union_gene_all('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_genes.txt')
         
interval = {}
for c in genes:
    interval[c] = []
    for gene_name in genes[c]:
        overlap = union_gene_all[union_gene_all['gene_name'] == gene_name][0]
        chro = overlap['chr'].lstrip('chr')
        start = overlap['start']
        end = overlap['end']
        interval[c].append((chro , start - 3000000 , end + 3000000 , gene_name))
        

CCS = np.loadtxt(os.path.join(PCFolder , 'CCS_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT5 = np.loadtxt(os.path.join(PCFolder , 'NT5_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
NT6 = np.loadtxt(os.path.join(PCFolder , 'NT6_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F35 = np.loadtxt(os.path.join(PCFolder , 'F35_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)
F40 = np.loadtxt(os.path.join(PCFolder , 'F40_Traditonal_PC_200K_Compartment_200K.txt') , dtype = sig_type)


size = (12, 10)
Left1 = 0.2 ; HB = 0.3 ; width1 = 0.45 ; HH = 0.3
Left2 = 0.66 ; width2 = 0.3 

Outfolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\'


for c in genes:
# for c in ['c3']:
    OutFil = c + '_Compartment_DEG_Plot.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    print (c , len(genes[c]))
    for i in interval[c]:
        g = i[0]
        start = i[1] // R
        end = i[2] // R
        name = c + '_' + i[3] + '_Chr' + g
        gene = union_gene_all[union_gene_all['gene_name'] == i[3]][0]
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
        ax = fig.add_axes([Left1, HB , width1 , 0])
        xtick = [start * R , i[1] + 3000000 , i[2] - 3000000 , end * R]
        labels = [properU(xtick[0]) , 'S' , 'E' , properU(xtick[-1])]
        ax.set_xticks(xtick)
        ax.set_xticklabels(labels)

        ax.set_xlabel(name , fontsize=20)
        ax.set_xlim((start * R , end * R ))
        ax0 = add_ax(sig5 , [Left1, HB , width1 , 0.075] , 'F40')
        caxis_S_horizontal(ax0)
        ax1 = add_ax(sig4 , [Left1, HB + 0.075 , width1 , 0.075] , 'F35')
        caxis_S_horizontal(ax1)
        ax2 = add_ax(sig3 , [Left1, HB + 0.15 , width1 , 0.075] , 'NT6')
        caxis_S_horizontal(ax2)
        ax3 = add_ax(sig2 , [Left1, HB + 0.225, width1 , 0.075] , 'NT5')
        caxis_S_horizontal(ax3)
        ax4 = add_ax(sig1 , [Left1, HB + 0.3, width1 , 0.075] , 'CCS')
        caxis_S_horizontal(ax4)
        
        v = np.log2(np.array([gene['CCS_1'] , gene['CCS_2'] , gene['CCS_3']]) + 1)
        w = np.log2(np.array([gene['NT5_1'] , gene['NT5_2'] , gene['NT5_3']]) + 1)
        x = np.log2(np.array([gene['NT6_1'] , gene['NT6_2'] , gene['NT6_3']]) + 1)
        y = np.log2(np.array([gene['F35_1'] , gene['F35_2'] , gene['F35_3']]) + 1)
        z = np.log2(np.array([gene['F40_1'] , gene['F40_2'] , gene['F40_3']]) + 1)
        ax5 = fig.add_axes([Left2, HB , width2 , HH + 0.075])
        ax5.scatter([1,1,1] , v , s = 700)
        ax5.scatter([2,2,2] , w , s = 700)
        ax5.scatter([3,3,3] , x , s = 700)
        ax5.scatter([4,4,4] , y , s = 700)
        ax5.scatter([5,5,5] , z , s = 700)
        ax5.scatter([1,2,3,4,5] , [v.mean() , w.mean() , x.mean() , y.mean() , z.mean()] , marker = '_' , s = 1000 , c = 'black')
        ax5.set_xlim((0.5,6))
        ax5.set_ylim((0,10))
        ax5.set_xticks([1,2,3,4,5])
        ax5.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40'] ,fontsize = 20)
#        ax5.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
        
            
        pp.savefig(fig)
        plt.close(fig)
    pp.close()

