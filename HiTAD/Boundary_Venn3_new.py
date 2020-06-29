# -*- coding: utf-8 -*-
"""
Created on Fri May 01 10:57:44 2020

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
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles



# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['skyblue' , 'k' , 'yellow'])
#my_cmap.set_bad('#2672a1')

IS_type = np.dtype({'names':['chr' , 'IS'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float]})
 
union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
         
        
        
        
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
    ax.set_ylabel('Insulation Score')
    ax.set_xlim((0.5 , 4.5))
    ax.set_ylim((0 , 5))
    return fig
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()



def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
 
def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
        for i in gene:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()      

    
ins_type = np.dtype({'names' : ['chr' , 'site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'] , 
                     'formats' : ['S8' , np.int , np.float , np.float , np.float , np.float]})
                                        
                    
CCS_NT_nofESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\CCS_NT_nofESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_NT_nofESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_NT_nofESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
CCS_noNT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\CCS_noNT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_noNT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_noNT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)

CCS_noNT_nofESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\CCS_noNT_nofESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_NT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_NT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
CCS_NT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\CCS_NT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)


CCS_NT_nofESC = []
for i in CCS_NT_nofESC_ins:
    if (i['CCS'] != 0) and (i['NT5'] / i['fESC'] > 1) and (i['NT6'] / i['fESC'] > 1) and (i['CCS'] / i['fESC'] > 1):
        CCS_NT_nofESC.append(i)
CCS_NT_nofESC = np.array(CCS_NT_nofESC , dtype = ins_type)

noCCS_NT_nofESC = []
for i in noCCS_NT_nofESC_ins:
    if (i['NT5'] != 0) and (i['NT6'] != 0) and (i['NT5'] / i['CCS'] > 1) and (i['NT6'] / i['CCS'] > 1) and (i['NT5'] / i['fESC'] > 1)  and (i['NT6'] / i['fESC'] > 1):
        noCCS_NT_nofESC.append(i)
noCCS_NT_nofESC = np.array(noCCS_NT_nofESC , dtype = ins_type)


CCS_noNT_fESC = []
for i in CCS_noNT_fESC_ins:
    if (i['CCS'] != 0) and (i['fESC'] != 0) and (i['fESC'] / i['NT5'] > 1) and (i['fESC'] / i['NT6'] > 1) and (i['CCS'] / i['NT5'] > 1) and (i['CCS'] / i['NT6'] > 1):
        CCS_noNT_fESC.append(i)
CCS_noNT_fESC = np.array(CCS_noNT_fESC , dtype = ins_type)
    
noCCS_noNT_fESC = []
for i in noCCS_noNT_fESC_ins:
    if (i['fESC'] != 0) and (i['fESC'] / i['NT5'] > 1) and (i['fESC'] / i['NT6'] > 1) and (i['fESC'] / i['CCS'] > 1):
        noCCS_noNT_fESC.append(i)
noCCS_noNT_fESC = np.array(noCCS_noNT_fESC , dtype = ins_type)



CCS_noNT_nofESC = []
for i in CCS_noNT_nofESC_ins:
    if (i['CCS'] != 0) and (i['CCS'] / i['NT5'] > 1) and (i['CCS'] / i['NT6'] > 1) and (i['CCS'] / i['fESC'] > 1):
        CCS_noNT_nofESC.append(i)
CCS_noNT_nofESC = np.array(CCS_noNT_nofESC , dtype = ins_type)

noCCS_NT_fESC = []
for i in noCCS_NT_fESC_ins:
    if (i['fESC'] != 0) and (i['NT5'] != 0) and (i['NT6'] != 0) and (i['NT5'] / i['CCS'] > 1) and (i['NT6'] / i['CCS'] > 1) and (i['fESC'] / i['CCS'] > 1):
        noCCS_NT_fESC.append(i)
noCCS_NT_fESC = np.array(noCCS_NT_fESC , dtype = ins_type)

CCS_NT_fESC = []
for i in CCS_NT_fESC_ins:
    if (i['CCS'] != 0) and (i['NT5'] != 0) and (i['NT6'] != 0) and (i['fESC'] != 0):
        CCS_NT_fESC.append(i)
CCS_NT_fESC = np.array(CCS_NT_fESC , dtype = ins_type)


Abc = len(CCS_noNT_nofESC)
aBc = len(noCCS_NT_nofESC)
ABc = len(CCS_NT_nofESC)
abC = len(noCCS_noNT_fESC)
AbC = len(CCS_noNT_fESC)
aBC = len(noCCS_NT_fESC)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC


fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('CCS', 'NT' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig_DY\\Plot\\Fig3_DY\\CCS_NT(common)fESC_Boundary_Venn3.pdf')


boundaries = {'CCS_noNT_nofESC' : CCS_noNT_nofESC , 'CCS_noNT_fESC' : CCS_noNT_fESC , 'CCS_NT_nofESC' : CCS_NT_nofESC,
              'noCCS_NT_fESC' : noCCS_NT_fESC , 'noCCS_NT_nofESC' : noCCS_NT_nofESC , 'noCCS_noNT_fESC' : noCCS_noNT_fESC ,
              'CCS_NT_fESC' : CCS_NT_fESC}



union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)



genes = {}

for cl in boundaries:
    genes[cl] = []
    for i in boundaries[cl]:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000 
        tmp_gene = union[union['chr'] == 'chr' + chro]        
        mask = (tmp_gene['gene_site'] >= start) & (tmp_gene['gene_site'] <= end)
        overlap = tmp_gene[mask]
        if overlap.size != 0 :
            for j in overlap:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
                
for k,v in genes.items():
    print k , len(v)
    Write2genes('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes_new\\' + k + '_genes.txt' , v)                
















fig1 = Box_plot([CCS_NT_nofESC['CCS'] , CCS_NT_nofESC['NT5'] , CCS_NT_nofESC['NT6'] , CCS_NT_nofESC['fESC']])
fig2 = Box_plot([noCCS_NT_nofESC['CCS'] , noCCS_NT_nofESC['NT5'] , noCCS_NT_nofESC['NT6'] , noCCS_NT_nofESC['fESC']])
fig3 = Box_plot([CCS_noNT_fESC['CCS'] , CCS_noNT_fESC['NT5'] , CCS_noNT_fESC['NT6'] , CCS_noNT_fESC['fESC']])
fig4 = Box_plot([noCCS_noNT_fESC['CCS'] , noCCS_noNT_fESC['NT5'] , noCCS_noNT_fESC['NT6'] , noCCS_noNT_fESC['fESC']])



fig1 = Box_plot([CCS_NT_nofESC_ins['CCS'] , CCS_NT_nofESC_ins['NT5'] , CCS_NT_nofESC_ins['NT6'] , CCS_NT_nofESC_ins['fESC']])
fig2 = Box_plot([noCCS_NT_nofESC_ins['CCS'] , noCCS_NT_nofESC_ins['NT5'] , noCCS_NT_nofESC_ins['NT6'] , noCCS_NT_nofESC_ins['fESC']])
fig3 = Box_plot([CCS_noNT_fESC_ins['CCS'] , CCS_noNT_fESC_ins['NT5'] , CCS_noNT_fESC_ins['NT6'] , CCS_noNT_fESC_ins['fESC']])
fig4 = Box_plot([noCCS_noNT_fESC_ins['CCS'] , noCCS_noNT_fESC_ins['NT5'] , noCCS_noNT_fESC_ins['NT6'] , noCCS_noNT_fESC_ins['fESC']])


run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig_DY\\Plot\\Fig3_DY\\CCS_NT_nofESC_Buondary_insulation_score.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig_DY\\Plot\\Fig3_DY\\noCCS_NT_nofESC_Buondary_insulation_score.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig_DY\\Plot\\Fig3_DY\\CCS_noNT_fESC_Buondary_insulation_score.pdf')
run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig_DY\\Plot\\Fig3_DY\\noCCS_noNT_fESC_Buondary_insulation_score.pdf')
