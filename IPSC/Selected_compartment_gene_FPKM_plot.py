# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 20:25:27 2021

@author: xxli
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


def get_union_gene_all_site_IPSC(Fil):
    union_type = ({'names':['gene_id' , 'gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3','F35_1','F35_2','F35_3','F40_1','F40_2','F40_3','MEF_1','MEF_2','IPSC_1','IPSC_2','E14_1','E14_2'],
                 'formats':['U64' , 'U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float , np.float, np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene


def Sort(a , s_1):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1]))
    return a

    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NTs':1 , 'fESC':2 , 'MEF':3 , 'IPSC':4 , 'E14':5}
cell = ['CCS','NTs','fESC','MEF','IPSC','E14']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
R = 200000
res = '200K'



sig_type = np.dtype({'names':['chr','score'],
                      'formats':['U8',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['U8',np.int , np.int , np.float]})
pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['U8' , np.float]})
gene_type = np.dtype({'names':['gene_names'] ,
                      'formats':['U64']})


Outfolder = 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S5\\'


c1 = 'Thsd7b Osgin2 Megf9 Hook1 Sema3d Hs3st1 Gabra4 Parm1 Arhgef5 Clec2h Zfp772 Hpgd Frk Caps2 Shisa6 Stxbp6 Galc Serpinb6b Serpinb1a Serpinb9b Serpinb9 Serpinb9e Serpinb6a Cpox Klhl28'
c2 = 'Bmpr1b Zar1 Zeb2 Ltbp1'
c3 = 'Sulf1 Eya1 Gm10561 Plcl1 Bard1 En1 Tagln2 Frmd4a Hoxd13 Dclk1 Rarres1 Plag1 Nfib Mfsd2a Fhl3 Laptm5 Megf6 Gm15772 Nsg1 Mn1 Cmklr1 Tmem119 Cit Mmp17 Flt1 Col1a2 Plxna4 Aqp1 Rad51ap1 Cacng7 Ehd2 Tshz3 Etfb Iglon5 Rgma Fam174b Dchs1 Syt17 Sfrp1 Chst2 Tgfbr2 Bicc1 Lgr5 Nipal4 Adam19 Cdc42ep4 Klhl29 Mirg Nid1 Id4 Cxcl14 Tgfbi Nid2 Gng2 Anxa8 Gzme Mcpt8 Ebf2 Enox1 Slc1a3 Ank Apol10b Prr5 Prickle1 Cebpd Hunk Erg Zfp521 Mamdc2 Sorbs1 Aldh18a1 Gsto1 Sec16b Tnn Snhg8 Gm5292 Cd38 Wbscr17 Pdzrn3 Lrig3 Nme2 Sox11 Lmo7 Fndc1 Sema6a Fam107b Itga8 Apbb1ip Fli1 Ust Tmem200a Gdf10'
c4 = 'Trim28 Peg3 Peg13 Snai2 Psg16 Psg19'
c5 = 'Col3a1 Tmeff2 Kcne4 Hmcn1 Nuf2 Rgs4 Rgs5 Jag1 Pcdh18 Glrb Col11a1 Arsj Tox Adamtsl1 Ldb2 Slit2 Gng11 Ptn Gas2 Fancf Gm15645 Nrg1 Tll1 Rnf150 Cdh11 Ets1 Lum Adra1b Osr1 Twist1 Flrt2 Gli3 Vcan Pde4d Plk2 Kctd12 Gpc6 Ghr Rspo2 Tnfrsf11b Enpp2 Abi3bp Epha3 Efna5 Slc8a1 Setbp1 Acta2 Stambpl1 Fas Ifit2 Ch25h Arhgap24 Man2a1 Cdon Shox2 Tbx15 Rbms3'
c6 = 'CAV2 CAV1 APOH THBS1 Trp63 DLX5 DLX6'
c7 = 'Nr5a2 Ildr2 F11r Irf6 Bend7 Slc33a1 Sfrp2 Ovgp1 Mob3b Kit Thap6 Idua Fgfrl1 Zfp113 N4bp2l1 C1galt1 Gimap9 Tuba8 Aicda Rimklb Tmem86a Stxbp2 Zfp709 Zfp617 Cyp4f18 Zfp612 Slc37a2 Gramd1b Slc35f2 Dmxl2 Klhl31 Foxl2os Foxl2 Pik3cb Plxnb1 Gpd1l Cmtm8 Ncoa7 Bves Slc17a8 Nr1h4 Tns3 Fam161a Usp43 Slc9a3r1 Synj2bp Dpf3 Slc12a7 Naip1 Ocln Ndrg2 Card6 Csf2rb Lnpep Slc37a1 Rsph1 H2-Ab1 Tap2 Slc25a23 Crb3 Rnf125 Tshz1 Zadh2 Pip5k1b Entpd1 Pyroxd2 Igsf11'
c8 = 'Sox17 Lin28b'

c1 = c1.split()
c2 = c2.split()
c3 = c3.split()
c4 = c4.split()
c5 = c5.split()
c6 = c6.split()
c7 = c7.split()
c8 = c8.split()

c6 = [x.capitalize() for x in c6]


##cluster_related_genes
c1_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_A_B_Donor_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c2_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_A_B_Method_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c3_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_B_A_Donor_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c4_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_B_A_Method_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c5_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_A_B_Donor_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c6_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_A_B_Method_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c7_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_B_A_Donor_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))
c8_genes = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_B_A_Method_specific_genes.txt' , skiprows=1 , dtype = gene_type , usecols=(0))





#interval = [('3' , 32000000 , 40000000)]
genes = {'NT_A_B_Donor':c1 ,'NT_A_B_Methods':c2 , 'NT_B_A_Donor':c3 , 'NT_B_A_Methods':c4 , 
         'IPSC_A_B_Donor':c5 , 'IPSC_A_B_Methods':c6 , 'IPSC_B_A_Donor':c7 , 'IPSC_B_A_Methods':c8}

related_genes = {'NT_A_B_Donor':c1_genes ,'NT_A_B_Methods':c2_genes , 'NT_B_A_Donor':c3_genes , 'NT_B_A_Methods':c4_genes , 
         'IPSC_A_B_Donor':c5_genes , 'IPSC_A_B_Methods':c6_genes , 'IPSC_B_A_Donor':c7_genes , 'IPSC_B_A_Methods':c8_genes}

union_gene_all = get_union_gene_all_site_IPSC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\gene_expression\\NT_IPSC_all_genes.txt')
         
interval = {}
for c in genes:
    interval[c] = []   
    if 'Donor' in c:
        for gene_name in genes[c]:
            if gene_name in related_genes[c]['gene_names']:
                overlap = union_gene_all[union_gene_all['gene_name'] == gene_name][0]
                chro = overlap['chr'].lstrip('chr')
                start = overlap['start']
                end = overlap['end']
                interval[c].append((chro , start - 3000000 , end + 3000000 , gene_name))
    elif 'Methods' in c:
        for gene_name in related_genes[c]['gene_names']:
            overlap = union_gene_all[union_gene_all['gene_name'] == gene_name][0]
            chro = overlap['chr'].lstrip('chr')
            start = overlap['start']
            end = overlap['end']
            interval[c].append((chro , start - 3000000 , end + 3000000 , gene_name))        
        
            

CCS = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
CCS_1 = Sort(list(CCS) , 0)
CCS = np.array(CCS_1 , dtype = CCS.dtype)
NTs = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NTs_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
NTs = NTs[NTs['chr'] != 'X']
NTs_1 = Sort(list(NTs) , 0)
NTs = np.array(NTs_1 , dtype = CCS.dtype)
fESC = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\fESC_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']
fESC_1 = Sort(list(fESC) , 0)
fESC = np.array(fESC_1 , dtype = CCS.dtype)

MEF = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
MEF_1 = Sort(list(MEF) , 0)
MEF = np.array(MEF_1 , dtype = CCS.dtype)
IPS_P3 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
IPS_P3_1 = Sort(list(IPS_P3) , 0)
IPS_P3 = np.array(IPS_P3_1 , dtype = CCS.dtype)
E14 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']
E14_1 = Sort(list(E14) , 0)
E14 = np.array(E14_1 , dtype = CCS.dtype)







#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 10)
Left1 = 0.2 ; HB = 0.2 ; width1 = 0.45 ; HH = 0.3
Left2 = 0.66 ; width2 = 0.3 




for c in genes:
    if 'IPSC_B_A_Methods' in c:
        OutFil = c + '_Compartment_gene_expression_1.pdf'
        pp = PdfPages(os.path.join(Outfolder , OutFil))
    # for c in ['c3']:
        print (c)
        for i in interval[c][:20]:
            g = i[0]
            start = i[1] // R
            end = i[2] // R
            name = c + '_' + i[3] + '_Chr' + g
            gene = union_gene_all[union_gene_all['gene_name'] == i[3]][0]
            tmp_CCS = CCS[CCS['chr'] == g]
            tmp_NTs = NTs[NTs['chr'] == g]
            tmp_fESC = fESC[fESC['chr'] == g]
            tmp_MEF = MEF[MEF['chr'] == g]
            tmp_IPSC = IPS_P3[IPS_P3['chr'] == g]
            tmp_E14 = E14[E14['chr'] == g]
            sig1 = tmp_CCS['pc'][start:end]
            sig2 = tmp_NTs['pc'][start:end]
            sig3 = tmp_fESC['pc'][start:end]
            sig4 = tmp_MEF['pc'][start:end]
            sig5 = tmp_IPSC['pc'][start:end]
            sig6 = tmp_E14['pc'][start:end]
            ##PCA Tracks
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left1, HB , width1 , 0])
            ax.set_xlim((start * R , end * R))
            xtick = [start * R , i[1] + 3000000 , i[2] - 3000000 , end * R]
            print (xtick)
            labels = [properU(xtick[0]) , 'S' , 'E' , properU(xtick[-1])]
            ax.set_xticks(xtick)            
            ax.set_xticklabels(labels)    
            ax.set_xlabel(name , fontsize=20)
            
            ax0 = add_ax(sig6 , [Left1, HB , width1 , 0.075] , 'E14')
            caxis_S_horizontal(ax0)
            ax1 = add_ax(sig5 , [Left1, HB + 0.075 , width1 , 0.075] , 'IPSC')
            caxis_S_horizontal(ax1)
            ax2 = add_ax(sig4 , [Left1, HB + 0.15 , width1 , 0.075] , 'MEF')
            caxis_S_horizontal(ax2)
            ax3 = add_ax(sig3 , [Left1, HB + 0.225, width1 , 0.075] , 'fESC')
            caxis_S_horizontal(ax3)
            ax4 = add_ax(sig2 , [Left1, HB + 0.3, width1 , 0.075] , 'NTs')
            caxis_S_horizontal(ax4)
            ax5 = add_ax(sig1 , [Left1, HB + 0.375, width1 , 0.075] , 'CCS')
            caxis_S_horizontal(ax5)
            
            u = np.log2(np.array([gene['CCS_1'] , gene['CCS_2'] , gene['CCS_3']]) + 1)
            v = np.log2(np.array([gene['NT5_1'] , gene['NT5_2'] , gene['NT5_3'] , gene['NT6_1'] , gene['NT6_2'] , gene['NT6_3']]) + 1)
            w = np.log2(np.array([gene['F35_1'] , gene['F35_2'] , gene['F35_3'] , gene['F40_1'] , gene['F40_2'] , gene['F40_3']]) + 1)
            x = np.log2(np.array([gene['MEF_1'] , gene['MEF_2']]) + 1)
            y = np.log2(np.array([gene['IPSC_1'] , gene['IPSC_2']]) + 1)
            z = np.log2(np.array([gene['E14_1'] , gene['E14_2']]) + 1)
            
            ax6 = fig.add_axes([Left2, HB , 0.3 , 0.4])
            ax6.scatter([1,1,1] , u , s = 700)
            ax6.scatter([2,2,2,2,2,2] , v , s = 700)
            ax6.scatter([3,3,3,3,3,3] , w , s = 700)
            ax6.scatter([4,4] , x , s = 700)
            ax6.scatter([5,5] , y , s = 700)
            ax6.scatter([6,6] , z , s = 700)
            
            ax6.scatter([1,2,3,4,5,6] , [u.mean() , v.mean() , w.mean() , x.mean() , y.mean() , z.mean()] , marker = '_' , s = 1000 , c = 'black')
            ax6.set_xlim((0.5,7))
            ax6.set_ylim((0,10))
            ax6.set_xticks([1,2,3,4,5,6])
            ax6.set_xticklabels(['CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC' , 'E14'] ,fontsize = 20)
    #        ax5.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
            
                
            pp.savefig(fig)
            plt.close(fig)
        pp.close()
            




