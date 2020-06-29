# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 20:46:25 2019

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

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


def get_union_gene(Fil):
    union_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                          'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union_gene = list(union_gene)
    union_gene = Sort(union_gene , 1 , 3)
    union_gene = np.array(union_gene , dtype = union_type)
    return union_gene
    
    
def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_name = i.strip().split('gene_name')[1].split('\"')[1]
            gene_id = i.strip().split('gene_id')[1].split('\"')[1]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf.close()
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf
    
    
def get_geneID(goFil):
    d_type = ({'names':['gene_id'],
               'formats':['S128']})
    goData = xlrd.open_workbook(goFil)
    sheet = goData.sheets()[0]
    pathway = sheet.col_values(1)
    p_value = sheet.col_values(4)
    goDatas = []
    for i in range(len(pathway)):
        if '~' not in pathway[i]:
            continue
        p_name = pathway[i].strip().split('~')[1]
        p_v = -np.log10(p_value[i])
        goDatas.append((p_name , p_v) )
    goDatas = np.array(goDatas , dtype = d_type)
    return goDatas

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a


#Get_genome    
gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})    
union_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                       'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

gtf_tmp = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
gtf = []
for i in gtf_tmp:
    gtf.append((i['gene_id'].split('.')[0] , i['gene_name'] , i['chr'] , i['strand'] , i['start'] , i['end']))
    
gtf = Sort(gtf , 2 , 4)
gtf = np.array(gtf , dtype = gtf_type)
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')  


Go_geneFil = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\Go_gene_ID.txt' , 'r')

for i in Go_geneFil:
    Go_CCS = i.strip().split(', ')
    break

for i in Go_geneFil:
    Go_ESC = i.strip().split(', ')
    break
Go_geneFil.close()
    
Go_gene_CCS = []
Go_gene_ESC = []    
for i in Go_CCS:
    gene_name = gtf[gtf['gene_id'] == i][0]['gene_name']
    gene = union_gene[union_gene['gene_name'] == gene_name][0]
    if (gene['CCS'] / gene['NT5'] > 1.5) and (gene['CCS'] / gene['NT6'] > 1.5) and (gene['CCS'] / gene['fESC'] > 1.5):
        Go_gene_CCS.append(gene)

Go_gene_CCS = np.array(Go_gene_CCS , dtype = union_type)

for i in Go_ESC:
    gene_name = gtf[gtf['gene_id'] == i][0]['gene_name']
    gene = union_gene[union_gene['gene_name'] == gene_name][0]
    if (gene['NT5'] / gene['CCS'] > 1.5) and (gene['NT6'] / gene['CCS'] > 1.5) and (gene['fESC'] / gene['CCS'] > 1.5):
        Go_gene_ESC.append(gene)
        
Go_gene_ESC = np.array(Go_gene_ESC , dtype = union_type)
        
Select_cluster2_gene = []
Select_cluster3_gene = []
Select_cluster2_genename = ['Ets1' , 'Runx1' , 'Runx2' , 'Irf2' , 'Tcf12' , 'Atf7']
Select_cluster3_genename = ['Sox2' , 'Sox3' , 'Pou5f1' , 'Klf5' , 'Ctcf' , 'Nanog']
Select_cluster2_genename.extend(list(Go_gene_CCS['gene_name']))
Select_cluster3_genename.extend(list(Go_gene_ESC['gene_name']))


for i in Select_cluster2_genename:
    gene = union_gene[union_gene['gene_name'] == i][0]
    Select_cluster2_gene.append(gene)
Select_cluster2_gene = np.array(Select_cluster2_gene , dtype = union_type)

for i in Select_cluster3_genename:
    gene = union_gene[union_gene['gene_name'] == i][0]
    Select_cluster3_gene.append(gene)
Select_cluster3_gene = np.array(Select_cluster3_gene , dtype = union_type)

w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\Select_cluster2_gene_expression_matrix.txt' , 'w')
for i in Select_cluster2_gene:
    w.writelines('\t'.join([str(np.log2(i['CCS'] + 1)) , str(np.log2(i['NT5'] + 1)) , str(np.log2(i['NT6'] + 1)) , str(np.log2(i['fESC'] + 1))]) + '\n')
w.close()

rownames(matrix) = c('Ets1','Runx1','Runx2','Irf2','Tcf12','Atf7','Rab20','Gm28635','Foxl2','Nceh1','Scd3','Ephx2','Samd8','Bmp6','Cyr61','Ntf3','Sema7a','Chst11','Ar','Ttll7')        
        
rownames(matrix) = c('Sox2','Sox3','Pou5f1','Klf5','Ctcf','Nanog','Ascl2','Pold2','Fgf10','Tipin','Foxm1','Sept9','Terf1','Ofd1','Utp14b','Usp44','Uchl5','Msh6','Spi1','Stn1')

        
        
        