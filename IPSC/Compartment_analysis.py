# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 19:52:19 2019

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

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['mediumblue' , 'white' , 'red'])
#my_cmap.set_bad('#2672a1')
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                      'formats':[np.float , np.float , np.float , np.float]})

locus_type = np.dtype({'names':['chr' , 'locus'] , 
                       'formats':['S8' , np.int]})

pc_NT_type = np.dtype({'names':['chr' , 'locus' , 'CCS' , 'NT5' , 'NT6' , 'fESC'] , 
                       'formats':['S8' , np.int , np.float , np.float , np.float , np.float]})
                    
pc_IPSC_type = np.dtype({'names':['chr' , 'locus' , 'MEF' , 'IPS_P3' , 'E14'] , 
                         'formats':['S8' , np.int , np.float , np.float , np.float]})
                    

def get_locus(A_B):
    locus_A_B = []
    for i in A_B:
        locus_A_B.append((i[-2] , i[-1]))
    locus_A_B = np.array(locus_A_B , dtype = locus_type)
    return locus_A_B
    
    
def get_Common_changed_PC(locus1 , locus2):
    Common = []
    loc1 = []
    loc2 = []
    for g in chroms:
        tmp1 = locus1[locus1['chr'] == g]
        tmp2 = locus2[locus2['chr'] == g]
        for i in tmp1:
            mask = (tmp2['locus'] == i['locus'])
            overlap = tmp2[mask]
            if overlap.size != 0:
                Common.append((i , overlap[0]))
            else:
                loc1.append((i['chr'] , i['locus']))
        for i in tmp2:
            mask = (tmp1['locus'] == i['locus'])
            overlap = tmp1[mask]
            if overlap.size == 0:
                loc2.append((i['chr'] , i['locus']))
            else:
                pass
    loc1 = np.array(loc1 , dtype = locus_type)
    loc2 = np.array(loc2 , dtype = locus_type)      
     
    return Common , loc1 , loc2

def get_union_changed_PC(locus1 , locus2):
    union = []
    for g in chroms:
        tmp1 = locus1[locus1['chr'] == g]
        tmp2 = locus2[locus2['chr'] == g]
        for i in tmp1:
            mask = (tmp2['locus'] == i['locus'])
            overlap = tmp2[mask]
            if overlap.size != 0:
                union.append((i , overlap[0]))
            else:
                union.append((i))
        for i in tmp2:
            mask = (tmp1['locus'] == i['locus'])
            overlap = tmp1[mask]
            if overlap.size != 0:
                pass
            else:
                union.append((i))
    return union

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

def get_compartment_genes(compartment , union):
    gene = []
    for i in compartment:
        g = i[0]
        start = i[1]
        end = i[1] + res
        tmp = union[union['chr'] == 'chr' + g]
        mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
        overlap = tmp[mask]
        if overlap.size > 0:
            for j in overlap:
                gene.append(j)
    gene = np.array(gene , dtype = union.dtype)            
    return gene
    
    
def Box_plot(data , compartment_numbers , gene_numbers):                
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
    ax.set_ylabel('compartment_numbers: ' + str(compartment_numbers) + ' , gene_numbers: ' + str(gene_numbers))
    ax.set_xlim((0.5 , 4.5))
    ax.set_ylim((-0.1 , 50))
    return fig



def IPSC_Box_plot(data , compartment_numbers , gene_numbers):                
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
    ax.set_ylabel('compartment_numbers: ' + str(compartment_numbers) + ' , gene_numbers: ' + str(gene_numbers))
    ax.set_xlim((0.5 , 3.5))
    ax.set_ylim((4 , 19))
    return fig


def get_FPKM_matrix(gene):
    matrix = np.zeros((len(gene) , 4))
    for i in range(len(gene)):
        matrix[i][0] = np.log2(gene[i]['CCS'] + 1)
        matrix[i][1] = np.log2(gene[i]['NT5'] + 1)
        matrix[i][2] = np.log2(gene[i]['NT6'] + 1)
        matrix[i][3] = np.log2(gene[i]['fESC'] + 1)
    return matrix
    
def get_FPKM3_matrix(gene):
    matrix = np.zeros((len(gene) , 3))
    for i in range(len(gene)):
        matrix[i][0] = np.log2(gene[i]['CCS'] + 1)
        matrix[i][1] = np.log2((gene[i]['NT5'] + gene[i]['NT6']) / 2 + 1)
        matrix[i][2] = np.log2(gene[i]['fESC'] + 1)
    return matrix
    
def get_IPSC_FPKM_matrix(gene):
    matrix = np.zeros((len(gene) , 3))
    for i in range(len(gene)):
        matrix[i][0] = np.log2(gene[i]['MEF'] + 1)
        matrix[i][1] = np.log2(gene[i]['IPS_P3'] + 1)
        matrix[i][2] = np.log2(gene[i]['E14'] + 1)
    return matrix
    
    
def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

def K_means_cluster(matrix_0,matrix,n , reverse):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1] , reverse=reverse))
    return matrix_1

    
def plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    vmax = np.percentile(matrix,95)
    vmin = np.percentile(matrix,5)
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']

    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('genes',fontsize = 30)
    plt.title('Zscore ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig
    
def IPSC_plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    vmax = np.percentile(matrix,95)
    vmin = np.percentile(matrix,5)
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['MEF' , 'IPS_P3' , 'E14']

    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('genes',fontsize = 30)
    plt.title('Zscore ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig
    
    
def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a
	
 
def Sort_array(matrix , col):
    '''
    a: list of needed to be sorted
    '''
    index = np.lexsort([matrix[:,col]])
    matrix = matrix[index]
    return matrix   
 
def rset_classify(matrix , order):
    matrix_new = matrix[matrix[:,-1] == order[0]]
    for i in order[1:]:
        m = matrix[matrix[:,-1] == i]
        matrix_new = np.array(np.vstack((matrix_new ,m)))

    return matrix_new
    
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



    
#Compartment                    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']


MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


A_B = []
B_A = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and ((tmp_NT5[i]['pc'] < 0 and tmp_NT6[i]['pc'] < 0) or (tmp_NT5[i]['pc'] > 0 and tmp_NT6[i]['pc'] > 0)) and (tmp_fESC[i]['pc'] < 0):
            A_B.append((g , int(i * res) , tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc']))
        elif (tmp_CCS[i]['pc'] < 0) and ((tmp_NT5[i]['pc'] < 0 and tmp_NT6[i]['pc'] < 0) or (tmp_NT5[i]['pc'] > 0 and tmp_NT6[i]['pc'] > 0)) and (tmp_fESC[i]['pc'] > 0):
            B_A.append((g , int(i * res) , tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc']))

A_B = np.array(A_B , dtype = pc_NT_type)
B_A = np.array(B_A , dtype = pc_NT_type)


##IPSC
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            IPSC_A_B.append((g , int(i * res) , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            IPSC_B_A.append((g , int(i * res) , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))

IPSC_A_B = np.array(IPSC_A_B , dtype = pc_IPSC_type)
IPSC_B_A = np.array(IPSC_B_A , dtype = pc_IPSC_type)


#common
common_A_B , nt_A_B , ipsc_A_B = get_Common_changed_PC(A_B , IPSC_A_B)
common_B_A , nt_B_A , ipsc_B_A = get_Common_changed_PC(B_A , IPSC_B_A)
    
Common_A_B = []
Common_B_A = []


for i in common_A_B:
    Common_A_B.append((i[0][0] , i[0][1]))

for i in common_B_A:
    Common_B_A.append((i[0][0] , i[0][1]))
    
Common_A_B = np.array(Common_A_B , dtype = locus_type)
Common_B_A = np.array(Common_B_A , dtype = locus_type)


data = {'NT_A_B':A_B , 'NT_B_A':B_A , 'NT_noIPSC_A_B':nt_A_B , 'NT_noIPSC_B_A':nt_B_A , 
        'IPSC_A_B':IPSC_A_B , 'IPSC_B_A':IPSC_B_A , 'noNT_IPSC_A_B':ipsc_A_B , 'noNT_IPSC_B_A':ipsc_B_A , 
        'NT_IPSC_common_A_B':Common_A_B , 'NT_IPSC_common_B_A':Common_B_A}

for k in data:
    print k , len(data[k])
    out = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\Fig3E_dynamic_genes_NT_vs_IPSC\\' + k + '.txt' , 'w')
    for i in data[k]:
        out.writelines(i['chr'] + '\t' + str(i['locus']) + '\n')
    out.close()


#common_resist        
resist_A_B = []
for i in common_A_B:
    if ((i[0][3] > 0) or (i[0][4] > 0)) and (i[1][3] > 0):
        resist_A_B.append((i[0][0] , i[0][1]))
        
resist_B_A = [] 
for i in common_B_A:
    if ((i[0][3] < 0) or (i[0][4] < 0)) and (i[1][3] < 0):
        resist_B_A.append((i[0][0] , i[0][1]))
 
Common_resist_A_B = np.array(resist_A_B , dtype = locus_type)
Common_resist_B_A = np.array(resist_B_A , dtype = locus_type)        
    
#respective_resist
NT_resist_A_B = []
NT_resist_B_A = []
IPSC_resist_A_B = []
IPSC_resist_B_A = []

for i in A_B:
    if (i['NT5'] > 0) and (i['NT6'] > 0):
        NT_resist_A_B.append(i)
for i in B_A:
    if (i['NT5'] < 0) and (i['NT6'] < 0):
        NT_resist_B_A.append(i)
        
for i in IPSC_A_B:
    if i['IPS_P3'] > 0:
        IPSC_resist_A_B.append(i)
for i in IPSC_B_A:
    if i['IPS_P3'] < 0:
        IPSC_resist_B_A.append(i)
        
NT_resist_A_B = np.array(NT_resist_A_B , dtype = A_B.dtype)
NT_resist_B_A = np.array(NT_resist_B_A , dtype = B_A.dtype)
IPSC_resist_A_B = np.array(IPSC_resist_A_B , dtype = IPSC_A_B.dtype)
IPSC_resist_B_A = np.array(IPSC_resist_B_A , dtype = IPSC_B_A.dtype)

common_resist_A_B , nt_resist_A_B , ipsc_resist_A_B = get_Common_changed_PC(NT_resist_A_B , IPSC_resist_A_B)
common_resist_B_A , nt_resist_B_A , ipsc_resist_B_A = get_Common_changed_PC(NT_resist_B_A , IPSC_resist_B_A)


data = {'NT_resist_A_B':NT_resist_A_B , 'NT_resist_B_A':NT_resist_B_A , 'NT_noIPSC_resist_A_B':nt_resist_A_B , 'NT_noIPSC_resist_B_A':nt_resist_B_A , 
        'IPSC_resist_A_B':IPSC_resist_A_B , 'IPSC_resist_B_A':IPSC_resist_A_B , 'noNT_IPSC_resist_A_B':ipsc_resist_A_B, 'noNT_IPSC_resist_B_A':ipsc_resist_B_A, 
        'NT_IPSC_common_resist_A_B':Common_resist_A_B , 'NT_IPSC_common_resist_B_A':Common_resist_B_A}

for k in data:
    print k , len(data[k])
    out = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Respective_resist_compartment\\' + k + '.txt' , 'w')
    for i in data[k]:
        out.writelines(i['chr'] + '\t' + str(i['locus']) + '\n')
    out.close()




#union        
union_A_B = get_union_changed_PC(A_B , IPSC_A_B)
union_B_A = get_union_changed_PC(B_A , IPSC_B_A)
    

NT_resist_A_B = [] ; IPSC_resist_A_B = []
for i in union_A_B:
    if len(i) == 2:
        if (i[0][3] > 0) and (i[0][4] > 0) and (i[1][3] < 0):
            NT_resist_A_B.append((i[0][0] , i[0][1]))
    elif len(i) == 6:
        if (i[3] > 0) and (i[4] > 0):
            NT_resist_A_B.append((i[0] , i[1]))
    else:
        pass
        
for i in union_A_B:
    if len(i) == 2:
        if (i[0][3] < 0) and (i[0][4] < 0) and (i[1][3] > 0):
            IPSC_resist_A_B.append((i[0][0] , i[0][1]))
    elif len(i) == 5:
        if (i[3] > 0):
            IPSC_resist_A_B.append((i[0] , i[1]))
    else:
        pass
    
NT_resist_B_A = [] ; IPSC_resist_B_A = []
for i in union_B_A:
    if len(i) == 2:
        if (i[0][3] < 0) and (i[0][4] < 0) and (i[1][3] > 0):
            NT_resist_B_A.append((i[0][0] , i[0][1]))
    elif len(i) == 6:
        if (i[3] < 0) and (i[4] < 0):
            NT_resist_B_A.append((i[0] , i[1]))
    else:
        pass
        
for i in union_B_A:
    if len(i) == 2:
        if (i[0][3] > 0) and (i[0][4] > 0) and (i[1][3] < 0):
            IPSC_resist_B_A.append((i[0][0] , i[0][1]))
    elif len(i) == 5:
        if (i[3] < 0):
            IPSC_resist_B_A.append((i[0] , i[1]))
    else:
        pass
        

##gene_expression
union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_type_2 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
        
gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union_gene_1 = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        union_gene_1.append((i['gene_name'] , i['chr'] , (i['start'] + i['end']) / 2 , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union_gene = np.array(union_gene_1 , dtype = union_type)

IPSC_union_gene = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt' , dtype = union_type_1 , skiprows = 1)

IPSC_union = []
for i in IPSC_union_gene:
    gene_site = int((i['start'] + i['end']) / 2)
    if (i['MEF'] > 0) or (i['IPS_P3'] > 0) or (i['E14'] > 0):
        IPSC_union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC'] , i['MEF'] , i['IPS_P3'] , i['E14']))
IPSC_union = np.array(IPSC_union , dtype = union_type_2)




#Compartment_gene_expression
##respective
NT_A_B_gene = get_compartment_genes(A_B , union_gene)
IPSC_A_B_gene = get_compartment_genes(IPSC_A_B , IPSC_union)
NT_B_A_gene = get_compartment_genes(B_A , union_gene)
IPSC_B_A_gene = get_compartment_genes(IPSC_B_A , IPSC_union)


o1 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\NT_noIPSC_A_B_genes.txt' , 'w')
o2 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\noNT_IPSC_A_B_genes.txt' , 'w')
o3 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\NT_IPSC_common_A_B_genes.txt' , 'w')
o4 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\NT_noIPSC_B_A_genes.txt' , 'w')
o5 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\noNT_IPSC_B_A_genes.txt' , 'w')
o6 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\dynamic_genes_NT_vs_IPSC\\NT_IPSC_common_B_A_genes.txt' , 'w')


o1.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
o2.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
o3.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
o4.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
o5.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
o6.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')


for i in NT_A_B_gene:
    if i['gene_name'] not in IPSC_A_B_gene['gene_name']:
        o1.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')
    else:
        o3.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')
      
for i in IPSC_A_B_gene:
    if i['gene_name'] not in NT_A_B_gene['gene_name']:
        o2.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')


n1 = 0        
for i in NT_B_A_gene:
    if i['gene_name'] not in IPSC_B_A_gene['gene_name']:
        o4.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')
    else:
        n1 += 1
        o6.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')


n2 = 0     
for i in IPSC_B_A_gene:
    if i['gene_name'] not in NT_B_A_gene['gene_name']:
        n2 += 1
        o5.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')
        
o1.close()
o2.close()
o3.close()
o4.close()
o5.close()
o6.close()
            
##common
common_A_B_gene_NT = get_compartment_genes(Common_A_B , union_gene)
common_B_A_gene_NT = get_compartment_genes(Common_B_A , union_gene)
common_A_B_gene_IPSC = get_compartment_genes(Common_A_B , IPSC_union)
common_B_A_gene_IPSC = get_compartment_genes(Common_B_A , IPSC_union)

##common_resist
resist_A_B_gene_NT = get_compartment_genes(resist_A_B , union_gene)
resist_B_A_gene_NT = get_compartment_genes(resist_B_A , union_gene)
resist_A_B_gene_IPSC = get_compartment_genes(resist_A_B , IPSC_union)
resist_B_A_gene_IPSC_1 = get_compartment_genes(resist_B_A , IPSC_union)

resist_B_A_gene_IPSC = []
for i in resist_B_A_gene_IPSC_1:
    if (i['MEF'] < i['IPS_P3']) and (i['IPS_P3'] < i['E14']):
        resist_B_A_gene_IPSC.append(i)
        
resist_B_A_gene_IPSC = np.array(resist_B_A_gene_IPSC , dtype = resist_B_A_gene_IPSC_1.dtype)

##respective_resist
NT_resist_A_B_gene_1 = get_compartment_genes(NT_resist_A_B , union_gene)
IPSC_resist_A_B_gene_1 = get_compartment_genes(IPSC_resist_A_B , IPSC_union)
NT_resist_B_A_gene_1 = get_compartment_genes(NT_resist_B_A , union_gene)
IPSC_resist_B_A_gene_1 = get_compartment_genes(IPSC_resist_B_A , IPSC_union)

NT_resist_A_B_gene = []
IPSC_resist_A_B_gene = []
NT_resist_B_A_gene = []
IPSC_resist_B_A_gene = []

for i in NT_resist_A_B_gene_1:
    if (i['NT5'] > i['fESC']) or (i['NT6'] > i['fESC']):
        NT_resist_A_B_gene.append(i)
        
for i in NT_resist_B_A_gene_1:
    if (i['NT5'] < i['fESC']) or (i['NT6'] < i['fESC']):
        NT_resist_B_A_gene.append(i)        
        
for i in IPSC_resist_A_B_gene_1:
    if (i['IPS_P3'] > i['E14']):
        IPSC_resist_A_B_gene.append(i)
        
for i in IPSC_resist_B_A_gene_1:
    if (i['IPS_P3'] < i['E14']):
        IPSC_resist_B_A_gene.append(i) 
        
NT_resist_A_B_gene = np.array(NT_resist_A_B_gene , dtype = NT_resist_A_B_gene_1.dtype)
IPSC_resist_A_B_gene = np.array(IPSC_resist_A_B_gene , dtype = IPSC_resist_A_B_gene_1.dtype)  
NT_resist_B_A_gene = np.array(NT_resist_B_A_gene , dtype = NT_resist_B_A_gene_1.dtype)
IPSC_resist_B_A_gene = np.array(IPSC_resist_B_A_gene , dtype = IPSC_resist_B_A_gene_1.dtype)          

AB = [] ; BA = [] ; NT_AB = [] ; NT_BA = [] ; IPSC_AB = [] ; IPSC_BA = []

for i in NT_resist_A_B_gene:
    overlap = IPSC_resist_A_B_gene[IPSC_resist_A_B_gene['gene_name'] == i['gene_name']]
    if overlap.size != 0:
        AB.append(overlap[0])
        
        
for i in NT_resist_B_A_gene:
    overlap = IPSC_resist_B_A_gene[IPSC_resist_B_A_gene['gene_name'] == i['gene_name']]
    if overlap.size != 0:
        BA.append(overlap[0])
        
        
        
AB = np.array(AB , dtype = union_type_2)
BA = np.array(BA , dtype = union_type_2)

o1 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\NT_IPSC_A_B_common_genes.txt' , 'w')
o2 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\NT_IPSC_B_A_common_genes.txt' , 'w')
o3 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\NT_A_B_genes.txt' , 'w')
o4 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\NT_B_A_genes.txt' , 'w')
o5 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\IPSC_A_B_genes.txt' , 'w')
o6 = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\compartment_resist\\IPSC_B_A_genes.txt' , 'w')

o1.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14']) + '\n')
o2.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14']) + '\n')
o3.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
o4.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
o5.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'MEF' , 'IPS_P3' , 'E14']) + '\n')
o6.writelines('\t'.join(['Gene_name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'MEF' , 'IPS_P3' , 'E14']) + '\n')

for i in AB:
    gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
    o1.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['CCS']) , str(i['NT5']) , str(i['NT6'])  , str(i['fESC']) , str(i['MEF']) , str(i['IPS_P3']) , str(i['E14'])]) + '\n')
o1.close()

for i in BA:
    gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
    o2.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['CCS']) , str(i['NT5']) , str(i['NT6'])  , str(i['fESC']) , str(i['MEF']) , str(i['IPS_P3']) , str(i['E14'])]) + '\n')
o2.close()

for i in NT_resist_A_B_gene:
    if i['gene_name'] not in AB['gene_name']:
        gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
        o3.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['CCS']) , str(i['NT5']) , str(i['NT6'])  , str(i['fESC'])]) + '\n')
o3.close()

for i in NT_resist_B_A_gene:
    if i['gene_name'] not in BA['gene_name']:
        gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
        o4.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['CCS']) , str(i['NT5']) , str(i['NT6'])  , str(i['fESC'])]) + '\n')
o4.close()        

for i in IPSC_resist_A_B_gene:
    if i['gene_name'] not in AB['gene_name']:
        gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
        o5.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['MEF']) , str(i['IPS_P3']) , str(i['E14'])]) + '\n')
o5.close()    


for i in IPSC_resist_B_A_gene:
    if i['gene_name'] not in BA['gene_name']:
        gene = union_gene[union_gene['gene_name'] == i['gene_name']][0]
        o6.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , str(gene['end']) , str(i['MEF']) , str(i['IPS_P3']) , str(i['E14'])]) + '\n')
o6.close()    



#BoxPlot
##respective
fig_NT_A_B = Box_plot([NT_A_B_gene['CCS'] , NT_A_B_gene['NT5'] , NT_A_B_gene['NT6'] , NT_A_B_gene['fESC']] , len(A_B) , len(NT_A_B_gene))
fig_IPSC_A_B = IPSC_Box_plot([IPSC_A_B_gene['MEF'] , IPSC_A_B_gene['IPS_P3'] , IPSC_A_B_gene['E14']], len(IPSC_A_B) , len(IPSC_A_B_gene))

fig_NT_B_A = Box_plot([NT_B_A_gene['CCS'] , NT_B_A_gene['NT5'] , NT_B_A_gene['NT6'] , NT_B_A_gene['fESC']], len(B_A) , len(NT_B_A_gene))
fig_IPSC_B_A = IPSC_Box_plot([IPSC_B_A_gene['MEF'] , IPSC_B_A_gene['IPS_P3'] , IPSC_B_A_gene['E14']], len(IPSC_B_A) , len(IPSC_B_A_gene))

run_Plot(fig_NT_A_B , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\D_NT_compartment_A_B_gene_FPKM.pdf' )
run_Plot(fig_IPSC_A_B , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\D_IPSC_compartment_A_B_gene_FPKM.pdf')
run_Plot(fig_NT_B_A , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\D_NT_compartment_B_A_gene_FPKM.pdf')
run_Plot(fig_IPSC_B_A , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\D_IPSC_compartment_B_A_gene_FPKM.pdf')




##Common
fig_common_A_B_NT = Box_plot([common_A_B_gene_NT['CCS'] , common_A_B_gene_NT['NT5'] , common_A_B_gene_NT['NT6'] , common_A_B_gene_NT['fESC']] , len(Common_A_B) , len(common_A_B_gene_NT))
fig_common_A_B_IPSC = IPSC_Box_plot([common_A_B_gene_IPSC['MEF'] , common_A_B_gene_IPSC['IPS_P3'] , common_A_B_gene_IPSC['E14']] , len(Common_A_B) , len(common_A_B_gene_IPSC))

fig_common_B_A_NT = Box_plot([common_B_A_gene_NT['CCS'] , common_B_A_gene_NT['NT5'] , common_B_A_gene_NT['NT6'] , common_B_A_gene_NT['fESC']] , len(Common_B_A) , len(common_B_A_gene_NT))
fig_common_B_A_IPSC = IPSC_Box_plot([common_B_A_gene_IPSC['MEF'] , common_B_A_gene_IPSC['IPS_P3'] , common_B_A_gene_IPSC['E14']] , len(Common_B_A) , len(common_B_A_gene_IPSC))

run_Plot(fig_common_A_B_NT , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\NT_compartment_common_A_B_gene_FPKM.pdf' )
run_Plot(fig_common_A_B_IPSC , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\IPSC_compartment_common_A_B_gene_FPKM.pdf')
run_Plot(fig_common_B_A_NT , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\NT_compartment_common_B_A_gene_FPKM.pdf')
run_Plot(fig_common_B_A_IPSC , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\fig3_figs\\IPSC_compartment_common_B_A_gene_FPKM.pdf')




##Common_common_resist
fig_resist_A_B_NT = Box_plot([resist_A_B_gene_NT['CCS'] , resist_A_B_gene_NT['NT5'] , resist_A_B_gene_NT['NT6'] , resist_A_B_gene_NT['fESC']] , len(resist_A_B) , len(resist_A_B_gene_NT))
fig_resist_A_B_IPSC = IPSC_Box_plot([resist_A_B_gene_IPSC['MEF'] , resist_A_B_gene_IPSC['IPS_P3'] , resist_A_B_gene_IPSC['E14']]  , len(resist_A_B) , len(resist_A_B_gene_IPSC))


fig_resist_B_A_NT = Box_plot([resist_B_A_gene_NT['CCS'] , resist_B_A_gene_NT['NT5'] , resist_B_A_gene_NT['NT6'] , resist_B_A_gene_NT['fESC']]  , len(resist_B_A) , len(resist_B_A_gene_NT))
fig_resist_B_A_IPSC = IPSC_Box_plot([resist_B_A_gene_IPSC['MEF'] , resist_B_A_gene_IPSC['IPS_P3'] , resist_B_A_gene_IPSC['E14']]  , len(resist_B_A) , len(resist_B_A_gene_IPSC))


run_Plot(fig_resist_A_B_NT , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_common_common_resist_A_B_gene_FPKM.pdf' )
run_Plot(fig_resist_A_B_IPSC , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_common_common_resist_A_B_gene_FPKM.pdf')
run_Plot(fig_resist_B_A_NT , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_common_common_resist_B_A_gene_FPKM.pdf')
run_Plot(fig_resist_B_A_IPSC , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_common_common_resist_B_A_gene_FPKM.pdf')



#respective_resist
fig_NT_resist_A_B = Box_plot([NT_resist_A_B_gene['CCS'] , NT_resist_A_B_gene['NT5'] , NT_resist_A_B_gene['NT6'] , NT_resist_A_B_gene['fESC']] , len(NT_resist_A_B) , len(NT_resist_A_B_gene))
fig_IPSC_resist_A_B = IPSC_Box_plot([IPSC_resist_A_B_gene['MEF'] , IPSC_resist_A_B_gene['IPS_P3'] , IPSC_resist_A_B_gene['E14']], len(IPSC_resist_A_B) , len(IPSC_resist_A_B_gene))

fig_NT_resist_B_A = Box_plot([NT_resist_B_A_gene['CCS'] , NT_resist_B_A_gene['NT5'] , NT_resist_B_A_gene['NT6'] , NT_resist_B_A_gene['fESC']], len(NT_resist_B_A) , len(NT_resist_B_A_gene))
fig_IPSC_resist_B_A = IPSC_Box_plot([IPSC_resist_B_A_gene['MEF'] , IPSC_resist_B_A_gene['IPS_P3'] , IPSC_resist_B_A_gene['E14']], len(IPSC_resist_B_A) , len(IPSC_resist_B_A_gene))

run_Plot(fig_NT_resist_A_B , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_respective_resist_A_B_gene_FPKM.pdf' )
run_Plot(fig_IPSC_resist_A_B , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_respective_resist_A_B_gene_FPKM.pdf')
run_Plot(fig_NT_resist_B_A , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_respective_resist_B_A_gene_FPKM.pdf')
run_Plot(fig_IPSC_resist_B_A , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_respective_resist_B_A_gene_FPKM.pdf')




##Common_common_resist
resist_A_B_matrix_NT = get_FPKM_matrix(resist_A_B_gene_NT)
#sort_matrix = Sort_array(resist_A_B_matrix_NT , 0)
#resist_A_B_matrix_NT = Z_score(resist_A_B_matrix_NT)
resist_A_B_matrix_NT = K_means_cluster(resist_A_B_matrix_NT , resist_A_B_matrix_NT , 5 , False)
#resist_A_B_matrix_NT = rset_classify(resist_A_B_matrix_NT , [0,3,4,2,1])
fig1 = plot_heatmap(resist_A_B_matrix_NT[:, :4])    


resist_A_B_matrix_IPSC = get_IPSC_FPKM_matrix(resist_A_B_gene_IPSC)
#sort_matrix = Sort_array(resist_A_B_matrix_IPSC , 0)
#resist_A_B_matrix_IPSC = Z_score(resist_A_B_matrix_IPSC)
resist_A_B_matrix_IPSC = K_means_cluster(resist_A_B_matrix_IPSC , resist_A_B_matrix_IPSC , 6 , False)
resist_A_B_matrix_IPSC = rset_classify(resist_A_B_matrix_IPSC , [2,4,0,1,3,5])
fig2 = IPSC_plot_heatmap(resist_A_B_matrix_IPSC[:, :3])  


resist_B_A_matrix_NT = get_FPKM_matrix(resist_B_A_gene_NT)
#sort_matrix = Sort_array(resist_B_A_matrix_NT , 0)
#resist_B_A_matrix_NT = Z_score(resist_B_A_matrix_NT)
resist_B_A_matrix_NT = K_means_cluster(resist_B_A_matrix_NT , resist_B_A_matrix_NT , 5 , False)
#resist_B_A_matrix_NT = rset_classify(resist_B_A_matrix_NT , [0,3,2,4,1])
fig3 = plot_heatmap(resist_B_A_matrix_NT[:,:4])    


resist_B_A_matrix_IPSC = get_IPSC_FPKM_matrix(resist_B_A_gene_IPSC)
#sort_matrix = Sort_array(resist_B_A_matrix_IPSC , 0)
#resist_B_A_matrix_IPSC = Z_score(resist_B_A_matrix_IPSC)
resist_B_A_matrix_IPSC = K_means_cluster(resist_B_A_matrix_IPSC , resist_B_A_matrix_IPSC , 5 , False)
resist_B_A_matrix_IPSC = rset_classify(resist_B_A_matrix_IPSC , [1,4,0,3,2])
fig4 = IPSC_plot_heatmap(resist_B_A_matrix_IPSC[:,:3])  

run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_common_common_resist_A_B_gene_FPKM_heatmap.pdf' )
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_common_common_resist_A_B_gene_FPKM_heatmap.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_compartment_common_common_resist_B_A_gene_FPKM_heatmap.pdf')
run_Plot(fig4 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\IPSC_compartment_common_common_resist_B_A_gene_FPKM_heatmap.pdf')



##respective_resist
resist_A_B_matrix_NT = get_FPKM_matrix(NT_resist_A_B_gene)
#sort_matrix = Sort_array(resist_A_B_matrix_NT , 0)
#resist_A_B_matrix_NT = Z_score(resist_A_B_matrix_NT)
resist_A_B_matrix_NT = K_means_cluster(resist_A_B_matrix_NT , resist_A_B_matrix_NT , 5 , False)
#resist_A_B_matrix_NT = rset_classify(resist_A_B_matrix_NT , [0,3,4,2,1])
fig1 = plot_heatmap(resist_A_B_matrix_NT[:, :4])    


resist_A_B_matrix_IPSC = get_IPSC_FPKM_matrix(IPSC_resist_A_B_gene)
#sort_matrix = Sort_array(resist_A_B_matrix_IPSC , 0)
#resist_A_B_matrix_IPSC = Z_score(resist_A_B_matrix_IPSC)
resist_A_B_matrix_IPSC = K_means_cluster(resist_A_B_matrix_IPSC , resist_A_B_matrix_IPSC , 6 , False)
#resist_A_B_matrix_IPSC = rset_classify(resist_A_B_matrix_IPSC , [2,4,0,1,3,5])
fig2 = IPSC_plot_heatmap(resist_A_B_matrix_IPSC[:, :3])  


resist_B_A_matrix_NT = get_FPKM_matrix(NT_resist_B_A_gene)
#sort_matrix = Sort_array(resist_B_A_matrix_NT , 0)
#resist_B_A_matrix_NT = Z_score(resist_B_A_matrix_NT)
resist_B_A_matrix_NT = K_means_cluster(resist_B_A_matrix_NT , resist_B_A_matrix_NT , 5 , False)
#resist_B_A_matrix_NT = rset_classify(resist_B_A_matrix_NT , [0,3,2,4,1])
fig3 = plot_heatmap(resist_B_A_matrix_NT[:,:4])    


resist_B_A_matrix_IPSC = get_IPSC_FPKM_matrix(IPSC_resist_B_A_gene)
#sort_matrix = Sort_array(resist_B_A_matrix_IPSC , 0)
#resist_B_A_matrix_IPSC = Z_score(resist_B_A_matrix_IPSC)
resist_B_A_matrix_IPSC = K_means_cluster(resist_B_A_matrix_IPSC , resist_B_A_matrix_IPSC , 5 , False)
#resist_B_A_matrix_IPSC = rset_classify(resist_B_A_matrix_IPSC , [1,4,0,3,2])
fig4 = IPSC_plot_heatmap(resist_B_A_matrix_IPSC[:,:3])  




#Venn_Plot_compartment
fig = plt.figure(figsize = (10, 10))
venn2(((len(A_B) - len(common_A_B)) / len(CCS) , (len(IPSC_A_B) - len(common_A_B)) / len(CCS) , len(common_A_B) / len(CCS)) , set_labels=('NT_A_B', 'IPSC_A_B'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_IPSC_A_B_common_Venn2.pdf')

fig = plt.figure(figsize = (10, 10))
venn2(((len(B_A) - len(common_B_A)) / len(CCS) , (len(IPSC_B_A) - len(common_B_A)) / len(CCS) , len(common_B_A) / len(CCS)) , set_labels=('NT_B_A', 'IPSC_B_A'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig3_figs\\NT_IPSC_B_A_common_Venn2.pdf')

##Venn_Plot_resist_compartment
fig = plt.figure(figsize = (10, 10))
venn2((np.round((len(NT_resist_A_B) - len(common_resist_A_B)) / len(CCS) , 7) , np.round((len(IPSC_resist_A_B) - len(common_resist_A_B)) / len(CCS) , 7) , np.round(len(common_resist_A_B) / len(CCS) , 7)) , set_labels=('NT_resist_A_B', 'IPSC_resist_A_B'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\D_NT_IPSC_compartment_resist_A_B_Venn2.pdf')

fig = plt.figure(figsize = (10, 10))
venn2((np.round((len(NT_resist_B_A) - len(common_resist_B_A)) / len(CCS) , 7) , np.round((len(IPSC_resist_B_A) - len(common_resist_B_A)) / len(CCS) , 7) , np.round(len(common_resist_B_A) / len(CCS) , 7)) , set_labels=('NT_resist_B_A', 'IPSC_resist_B_A'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\D_NT_IPSC_compartment_resist_B_A_Venn2.pdf')


#Venn_Plot_gene_expression
fig = plt.figure(figsize = (10, 10))
venn2((109-15 , 163-15 , 15) , set_labels=('NT_resist_A_B', 'IPSC_resist_A_B'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_IPSC_resist_A_B_common_Venn2.pdf')

fig = plt.figure(figsize = (10, 10))
venn2((43-2 , 186-2 , 2) , set_labels=('NT_resist_B_A', 'IPSC_resist_B_A'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_IPSC_resist_B_A_common_Venn2.pdf')


o = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\compartment_resist_gene\\NT_Compartment_A_B_respective_resist_gene_names.txt' , 'w')
o.writelines('\t'.join(['Gene_names' , 'Gene_ID' , 'CCS_FPKM' , 'NT5_FPKM' , 'NT6_FPKM','fESC_FPKM']) + '\n')

for i in NT_resist_A_B_gene:
    gene_name = i['gene_name']
    gene_ID = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
    o.writelines('\t'.join([gene_name , gene_ID , str(i['CCS']) , str(i['NT5']) , str(i['NT6']) , str(i['fESC'])]) + '\n')
    
o.close()


