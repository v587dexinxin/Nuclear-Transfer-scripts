# -*- coding: utf-8 -*-
"""
Created on Thu Jan 09 15:43:30 2020

@author: han-luo
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:33:29 2019

@author: han-luo
"""
from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import copy
import scipy
from scipy import stats
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})

union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})


all_type = ({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_1' , 'CCS_2' , 'CCS_3' , 'NT5_1' , 'NT5_2' , 'NT5_3' , 'NT5_4' , 'NT6_1' , 'NT6_2' , 'NT6_3' , 'fESC_1' , 'fESC_2' , 'fESC_3'],
                 'formats':['S64' , 'S64', 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
                 
                 
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
    
def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf


 
def Box_plot_0(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'seagreen','linewidth':2},
            medianprops={'color':'seagreen','linewidth':2},
            capprops={'color':'seagreen','linewidth':2},
            whiskerprops={'color':'seagreen','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'chocolate','linewidth':2},
            medianprops={'color':'chocolate','linewidth':2},
            capprops={'color':'chocolate','linewidth':2},
            whiskerprops={'color':'chocolate','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7, 
            boxprops={'color': 'slateblue','linewidth':2},
            medianprops={'color':'slateblue','linewidth':2},
            capprops={'color':'slateblue','linewidth':2},
            whiskerprops={'color':'slateblue','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'deeppink','linewidth':2},
            medianprops={'color':'deeppink','linewidth':2},
            capprops={'color':'deeppink','linewidth':2},
            whiskerprops={'color':'deeppink','linewidth':2, 'linestyle':'--'})
#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[0] , data[3])[1]
    d4 = scipy.stats.ranksums(data[1] , data[2])[1]
    d5 = scipy.stats.ranksums(data[1] , data[3])[1]
    d6 = scipy.stats.ranksums(data[2] , data[3])[1]
    
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F40' ] , fontsize = 28)
    ax.set_xlabel('CCS_NT5:' + str(d1) + ',CCS_NT6:' + str(d2) + ',CCS_fESC:' + str(d3) + '\n' + 'NT5_NT6:' + str(d4) + ',NT5_fESC:' + str(d5) + ',NT6_fESC:' + str(d6))
    ax.set_xlim((0.5 , 4.5))
    ax.set_ylim((0 , 50))
    return fig


def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot([data[0]['CCS'] , data[1]['CCS'] , data[2]['CCS'] , data[3]['CCS']] , 
            positions=[1 , 6 , 11 , 16] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'seagreen','linewidth':1},
            medianprops={'color':'seagreen','linewidth':1},
            capprops={'color':'seagreen','linewidth':1},
            whiskerprops={'color':'seagreen','linewidth':1, 'linestyle':'--'})
    ax.boxplot([data[0]['NT5'] , data[1]['NT5'] , data[2]['NT5'] , data[3]['NT5']] ,
            positions=[2 , 7 , 12 , 17] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'chocolate','linewidth':1},
            medianprops={'color':'chocolate','linewidth':1},
            capprops={'color':'chocolate','linewidth':1},
            whiskerprops={'color':'chocolate','linewidth':1, 'linestyle':'--'})
    ax.boxplot([data[0]['NT6'] , data[1]['NT6'] , data[2]['NT6'] , data[3]['NT6']] ,
            positions=[3 , 8 , 13 , 18] , showfliers=False, widths = 0.7, 
            boxprops={'color': 'slateblue','linewidth':1},
            medianprops={'color':'slateblue','linewidth':1},
            capprops={'color':'slateblue','linewidth':1},
            whiskerprops={'color':'slateblue','linewidth':1, 'linestyle':'--'})
    ax.boxplot([data[0]['fESC'] , data[1]['fESC'] , data[2]['fESC'] , data[3]['fESC']] , 
            positions=[4 , 9 , 14 , 19] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'deeppink','linewidth':1},
            medianprops={'color':'deeppink','linewidth':1},
            capprops={'color':'deeppink','linewidth':1},
            whiskerprops={'color':'deeppink','linewidth':1, 'linestyle':'--'})
    ax.plot([5 , 5],[0 , 50], lw = 1.5, ls = '--', color = 'gray')
    ax.plot([10 , 10],[0 , 50], lw = 1.5, ls = '--', color = 'gray')
    ax.plot([15 , 15],[0 , 50], lw = 1.5, ls = '--', color = 'gray')
#    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
#    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
#    d3 = scipy.stats.ranksums(data[0] , data[3])[1]
#    d4 = scipy.stats.ranksums(data[1] , data[2])[1]
#    d5 = scipy.stats.ranksums(data[1] , data[3])[1]
#    d6 = scipy.stats.ranksums(data[2] , data[3])[1]
    
    ax.set_xticks([1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F40'  , '' , 'CCs' , 'NT5' , 'NT6' , 'F40'  , '' , 'CCs' , 'NT5' , 'NT6' , 'F40'  , '' , 'CCs' , 'NT5' , 'NT6' , 'F40'] , fontsize = 10)
    ax.set_xlabel('Repr , Partial , Over , Resist')
    ax.set_xlim((0 , 20))
    ax.set_ylim((0 , 50))
    return fig


def Write_comp_genename_ID(gene_cluster , gene_name , gene_ID):
    out1 = open(gene_name , 'w')
    out2 = open(gene_ID , 'w')
    for i in gene_cluster:
        gene_name = i['gene_name']
        out1.writelines(gene_name +'\n')
        gene_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
        out2.writelines(gene_id + '\n')
    out1.close()
    out2.close()

def Write_comp_genename_ID_expression(gene_cluster , outFil):
    out = open(outFil , 'w')
    out.writelines('\t'.join(['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT5_fPKM' , 'NT6_FPKM' , 'fESC_FPKM']) + '\n')
    for i in gene_cluster:
        gene_name = i['gene_name']
        gene_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
        gene = union_gene[union_gene['gene_name'] == gene_name][0]
        out.writelines(gene_id + '\t' + '\t'.join([str(a) for a in gene]) + '\n')
    out.close()


def Write_comp_genename_ID_expression_allreps(gene_cluster , outFil):
    out = open(outFil , 'w')
    out.writelines('\t'.join(['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_R1' , 'CCS_R2' , 'CCS_R3' , 'NT5_R1' , 'NT5_R2' , 'NT5_R3' , 'NT5_R4' , 'NT6_R1' , 'NT6_R2' , 'NT6_R3' , 'fESC_R1' , 'fESC_R2' , 'fESC_R3']) + '\n')
    for i in gene_cluster:
        gene_name = i['gene_name']
        gene = all_gene[all_gene['gene_name'] == gene_name][0]
        out.writelines('\t'.join([gene['gene_id'] , gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(gene['CCS_1']) , str(gene['CCS_2']) , str(gene['CCS_3']) , str(gene['NT5_1']) , str(gene['NT5_2']) , str(gene['NT5_3']) , str(gene['NT5_4']) , str(gene['NT6_1']) , str(gene['NT6_2']) , str(gene['NT6_3']) , str(gene['fESC_1']) , str(gene['fESC_2']) , str(gene['fESC_3'])]) + '\n')
    out.close()
    
    
    
def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
       
                      
c1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c7 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c8 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)


gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
all_gene = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt' , dtype = all_type , skiprows = 1)


union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)


c1_gene = [] ; c2_gene = [] ; c3_gene = []
c4_gene = [] ; c5_gene = [] ; c6_gene = []
c7_gene = [] ; c8_gene = []

for i in c1:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c1_gene.append(j)
c1_gene = np.array(c1_gene , dtype = union_type)     
       
for i in c2:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c2_gene.append(j)
c2_gene = np.array(c2_gene , dtype = union_type)         

for i in c3:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c3_gene.append(j)
c3_gene = np.array(c3_gene , dtype = union_type)         

for i in c4:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c4_gene.append(j)
c4_gene = np.array(c4_gene , dtype = union_type)      

for i in c5:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c5_gene.append(j)
c5_gene = np.array(c5_gene , dtype = union_type)     

for i in c6:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c6_gene.append(j)
c6_gene = np.array(c6_gene , dtype = union_type)     

         
for i in c7:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c7_gene.append(j)
c7_gene = np.array(c7_gene , dtype = union_type)     

for i in c8:
    chro = i['chr']
    start = i['start']
    end = i['end']
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c8_gene.append(j)
c8_gene = np.array(c8_gene , dtype = union_type)     

               

#fig1 = Box_plot([c1_gene['CCS'] , c1_gene['NT5'] , c1_gene['NT6'] , c1_gene['fESC']])
#fig2 = Box_plot([c2_gene['CCS'] , c2_gene['NT5'] , c2_gene['NT6'] , c2_gene['fESC']])
#fig3 = Box_plot([c3_gene['CCS'] , c3_gene['NT5'] , c3_gene['NT6'] , c3_gene['fESC']])
#fig4 = Box_plot([c4_gene['CCS'] , c4_gene['NT5'] , c4_gene['NT6'] , c4_gene['fESC']])
#fig5 = Box_plot([c5_gene['CCS'] , c5_gene['NT5'] , c5_gene['NT6'] , c5_gene['fESC']])
#fig6 = Box_plot([c6_gene['CCS'] , c6_gene['NT5'] , c6_gene['NT6'] , c6_gene['fESC']])
#fig7 = Box_plot([c7_gene['CCS'] , c7_gene['NT5'] , c7_gene['NT6'] , c7_gene['fESC']])
#fig8 = Box_plot([c8_gene['CCS'] , c8_gene['NT5'] , c8_gene['NT6'] , c8_gene['fESC']])
#
#
#run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster1.pdf')
#run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster2.pdf')
#run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster3.pdf')
#run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster4.pdf')
#run_Plot(fig5 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster5.pdf')
#run_Plot(fig6 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster6.pdf')
#run_Plot(fig7 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster7.pdf')
#run_Plot(fig8 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_cluster8.pdf')

fig1 = Box_plot([c1_gene , c2_gene , c3_gene , c4_gene])
fig2 = Box_plot([c5_gene , c6_gene , c7_gene , c8_gene])
run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_A-B.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_byself_new_gene_FPKM_boxplot_B-A.pdf')





Write_comp_genename_ID(c1_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_geneID.txt')


Write_comp_genename_ID(c2_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_geneID.txt')


Write_comp_genename_ID(c3_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_geneID.txt')


Write_comp_genename_ID(c4_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_geneID.txt')


Write_comp_genename_ID(c5_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_geneID.txt')


Write_comp_genename_ID(c6_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_geneID.txt')


Write_comp_genename_ID(c7_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_geneID.txt')


Write_comp_genename_ID(c8_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_geneName.txt'  ,
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_geneID.txt')


Write_comp_genename_ID_expression(c1_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_Gene.txt')

Write_comp_genename_ID_expression(c2_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_Gene.txt')

Write_comp_genename_ID_expression(c3_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_Gene.txt')

Write_comp_genename_ID_expression(c4_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_Gene.txt')

Write_comp_genename_ID_expression(c5_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_Gene.txt')

Write_comp_genename_ID_expression(c6_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_Gene.txt')

Write_comp_genename_ID_expression(c7_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_Gene.txt')

Write_comp_genename_ID_expression(c8_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_Gene.txt')



Write_comp_genename_ID_expression_allreps(c1_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c2_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c3_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c4_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c5_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c6_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c7_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_Gene_allreps.txt')

Write_comp_genename_ID_expression_allreps(c8_gene , 
'H:\\Workspace_New\\data\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_Gene_allreps.txt')






#-------------------------plot_compartment_classify_gene_number-------------------------------               
c1 = len(c1_gene) / 287 ; c2 = len(c2_gene) / 287 ; c3 = len(c3_gene) / 287  
c4 = len(c4_gene) / 287 ; c5 = len(c5_gene) / 1982 ; c6 = len(c6_gene) / 1982 
c7 = len(c7_gene) / 1982 ; c8 = len(c8_gene) / 1982

left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (6, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))

ax.bar(1 , c1 + c2 + c3 + c4 , color='green')
ax.bar(1 , c1 + c2 + c3 , color='red')
ax.bar(1 , c1 + c2 , color='orange')
ax.bar(1 , c1 , color='mediumblue')

ax.bar(2 , c5 + c6 + c7 + c8 , color='green')
ax.bar(2 , c5 + c6 + c7 , color='red')
ax.bar(2 , c5 + c6 , color='orange')
ax.bar(2 , c5 , color='mediumblue')



x = ['A-B','B-A']
ax.set_xticks([1 , 2])
ax.set_xticklabels(x,fontsize = 10)
ax.set_ylabel('Percent(%)',fontsize = 10)
ax.set_ylim((0.5 , 2.5))
ax.set_ylim((0 , 1))

        
plt.title('Compartment classified gene numbers',fontsize = 10)

    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_classified_gene_numbers.pdf') 


