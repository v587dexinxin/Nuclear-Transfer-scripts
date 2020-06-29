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


 
def Box_plot(data):                
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
    ax.set_ylim((0 , 40))
    return fig


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
       
                      
c1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)
c4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1 , 2) , dtype = data_type)

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)


c1_gene = []
c2_gene = []
c3_gene = []
c4_gene = []

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

fig1 = Box_plot([c1_gene['CCS'] , c1_gene['NT5'] , c1_gene['NT6'] , c1_gene['fESC']])
fig2 = Box_plot([c2_gene['CCS'] , c2_gene['NT5'] , c2_gene['NT6'] , c2_gene['fESC']])
fig3 = Box_plot([c3_gene['CCS'] , c3_gene['NT5'] , c3_gene['NT6'] , c3_gene['fESC']])
fig4 = Box_plot([c4_gene['CCS'] , c4_gene['NT5'] , c4_gene['NT6'] , c4_gene['fESC']])

run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_gene_FPKM_boxplot_cluster1.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_gene_FPKM_boxplot_cluster2.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_gene_FPKM_boxplot_cluster3.pdf')
run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\E_compartment_classify_gene_FPKM_boxplot_cluster4.pdf')

out1 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster4_geneName.txt' , 'w')
out2 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify\\Compartment_cluster4_geneID.txt' , 'w')
for i in c4_gene:
    gene_name = i['gene_name']
    out1.writelines(gene_name +'\n')
    gene_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
    out2.writelines(gene_id + '\n')
out1.close()
out2.close()