# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 11:16:49 2019

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



def IPSC_Box_plot(data):                
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

#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[1] , data[2])[1]

    
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels(['MEF' , 'IPS_P3' , 'E14'] , fontsize = 28)
    ax.set_xlabel('MEF_IPS_P3:' + str(d1) + ',MEF_E14:' + str(d2) + ',IPS_P3_E14:' + str(d3))
    ax.set_xlim((0.5 , 3.5))
#    ax.set_ylim((0 , 40))
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
       
pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})
union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
loc_type = np.dtype({'names':['chr' , 'start'] , 
                    'formats':['S8' , np.int]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000
                      
#Nuclear transfer                  
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']

#IPSC
MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_Compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_Compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_Compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']


#Nuclear transfer
A_B = []
B_A = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_NT5[i]['pc'] > 0) and (tmp_NT6[i]['pc'] > 0) and (tmp_fESC[i]['pc'] < 0):
            A_B.append((g , i))
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_NT5[i]['pc'] < 0) and (tmp_NT6[i]['pc'] < 0) and (tmp_fESC[i]['pc'] > 0):
            B_A.append((g , i))
A_B = np.array(A_B , dtype = loc_type)
B_A = np.array(B_A , dtype = loc_type)


##IPSC
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_IPS_P3[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            IPSC_A_B.append((g , i))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_IPS_P3[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            IPSC_B_A.append((g , i))

IPSC_A_B = np.array(IPSC_A_B , dtype = loc_type)
IPSC_B_A = np.array(IPSC_B_A , dtype = loc_type)


diff_A_B = [] ; diff_B_A =[] 
for g in chroms:
    tmp_A_B = A_B[A_B['chr'] == g]
    tmp_B_A = B_A[B_A['chr'] == g]
    tmp_IPSC_A_B = IPSC_A_B[IPSC_A_B['chr'] == g]
    tmp_IPSC_B_A = IPSC_B_A[IPSC_B_A['chr'] == g]
    for i in tmp_IPSC_A_B:
        if i['start'] not in tmp_A_B['start']:
            diff_A_B.append((g , i['start']))
        else:
            pass
    for i in tmp_IPSC_B_A:
        if i['start'] not in tmp_B_A['start']:
            diff_B_A.append((g , i['start']))
        else:
            pass
    
    
union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_type_2 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
              
       

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)

IPSC_union_gene = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt' , dtype = union_type_1 , skiprows = 1)

IPSC_union = []
for i in IPSC_union_gene:
    gene_site = int((i['start'] + i['end']) / 2)
    IPSC_union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , i['MEF'] , i['IPS_P3'] , i['E14']))
IPSC_union = np.array(IPSC_union , dtype = union_type_2)


c1 = [] ; c2 = []
for i in diff_A_B:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = IPSC_union[IPSC_union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c1.append(j)
c1 = np.array(c1 , dtype = union_type_2)     
       
for i in diff_B_A:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = IPSC_union[IPSC_union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c2.append(j)
c2 = np.array(c2 , dtype = union_type_2)

NT_c1 = [] ; NT_c2 = []
for i in diff_A_B:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c1.append(j)
NT_c1 = np.array(NT_c1 , dtype = union_type_2)  

for i in diff_B_A:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c2.append(j)
NT_c2 = np.array(NT_c2 , dtype = union_type_2)                 

                  

fig1 = Box_plot([NT_c1['CCS'] , NT_c1['NT5'] , NT_c1['NT6'] , NT_c1['fESC']])
fig2 = Box_plot([NT_c2['CCS'] , NT_c2['NT5'] , NT_c2['NT6'] , NT_c2['fESC']])

fig3 = IPSC_Box_plot([c1['MEF'] , c1['IPS_P3'] , c1['E14']])
fig4 = IPSC_Box_plot([c2['MEF'] , c2['IPS_P3'] , c2['E14']])


run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_ESC_noconsist_IPSCdiff_compartment_classify_gene_FPKM_boxplot_A_B.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_ESC_noconsist_IPSCdiff_compartment_classify_gene_FPKM_boxplot_B_A.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_ESC_noconsist_IPSCdiff_compartment_classify_gene_FPKM_boxplot_A_B.pdf')
run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_ESC_noconsist_IPSCdiff_compartment_classify_gene_FPKM_boxplot_B_A.pdf')


##consist or no consist

#Nuclear transfer
A_B = []
B_A = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_NT5[i]['pc'] < 0) and (tmp_NT6[i]['pc'] < 0) and (tmp_fESC[i]['pc'] < 0):
            A_B.append((g , i))
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_NT5[i]['pc'] > 0) and (tmp_NT6[i]['pc'] > 0) and (tmp_fESC[i]['pc'] > 0):
            B_A.append((g , i))
            

A_B = np.array(A_B , dtype = loc_type)
B_A = np.array(B_A , dtype = loc_type)


##IPSC
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_IPS_P3[i]['pc'] < 0) and (tmp_E14[i]['pc'] < 0):
            IPSC_A_B.append((g , i))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_IPS_P3[i]['pc'] > 0) and (tmp_E14[i]['pc'] > 0):
            IPSC_B_A.append((g , i))

IPSC_A_B = np.array(IPSC_A_B , dtype = loc_type)
IPSC_B_A = np.array(IPSC_B_A , dtype = loc_type)



    
union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_type_2 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
              
       

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)

IPSC_union_gene = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt' , dtype = union_type_1 , skiprows = 1)

IPSC_union = []
for i in IPSC_union_gene:
    gene_site = int((i['start'] + i['end']) / 2)
    IPSC_union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , i['MEF'] , i['IPS_P3'] , i['E14']))
IPSC_union = np.array(IPSC_union , dtype = union_type_2)




NT_c1 = [] ; NT_c2 = []
for i in A_B:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c1.append(j)
NT_c1 = np.array(NT_c1 , dtype = union_type_2)  

for i in B_A:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c2.append(j)
NT_c2 = np.array(NT_c2 , dtype = union_type_2)                 



c1 = [] ; c2 = []
for i in IPSC_A_B:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = IPSC_union[IPSC_union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c1.append(j)
c1 = np.array(c1 , dtype = union_type_2)     
       
for i in IPSC_B_A:
    chro = i[0]
    start = i[1] * 200000
    end = (i[1] + 1)* 200000
    tmp = IPSC_union[IPSC_union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c2.append(j)
c2 = np.array(c2 , dtype = union_type_2)




fig1 = Box_plot([NT_c1['CCS'] , NT_c1['NT5'] , NT_c1['NT6'] , NT_c1['fESC']])
fig2 = Box_plot([NT_c2['CCS'] , NT_c2['NT5'] , NT_c2['NT6'] , NT_c2['fESC']])

fig3 = IPSC_Box_plot([c1['MEF'] , c1['IPS_P3'] , c1['E14']])
fig4 = IPSC_Box_plot([c2['MEF'] , c2['IPS_P3'] , c2['E14']])

run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_ESC_consist_compartment_classify_gene_FPKM_boxplot_A_B.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_ESC_consist_compartment_classify_gene_FPKM_boxplot_B_A.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_ESC_consist_compartment_classify_gene_FPKM_boxplot_A_B.pdf')
run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\IPSC_ESC_consist_compartment_classify_gene_FPKM_boxplot_B_A.pdf')

out1 = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_compartment_A_B_geneNames.txt' , 'w')
out2 = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_compartment_A_B_geneID.txt' , 'w')

for i in c1:
    gene_name = i['gene_name']
    gene_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
    out1.writelines(gene_name + '\n')
    out2.writelines(gene_id + '\n')
out1.close()
out2.close()

