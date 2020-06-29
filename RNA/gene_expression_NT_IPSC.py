# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 21:46:32 2019

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


def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['gene_name' , 'chr' , 'strand' ,	'start' , 'end' , 'CCS_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'fESC_FPKM'	, 'MEF_normalized_signal' , 'IPS_P3_normalized_signal' , 'E14_normalized_signal']) + '\n')
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
    
    

data_type = np.dtype({'names':['geneID' , 'gene_name'] , 
                    'formats':['S64' , 'S64']})
                    
data_type_1 = np.dtype({'names':['geneID' , 'FPKM'] , 
                        'formats':['S64' , np.float]})

union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
                    
                    
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

MEF_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026297_MEF_1_1100252800520270.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))
MEF_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026298_MEF_2_2400252800520270.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))
IPS_P3_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026305_MEF_iPS_p3_1_1300252800519286.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))
IPS_P3_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026306_MEF_iPS_p3_2_2100252800519287.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))
E14_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026317_ES_1_1400252800520270.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))
E14_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\GSM2026318_ES_2_2200252800520270.txt' , dtype = data_type , skiprows = 13 , usecols = (11 , 12))

gene_id = []
for i in MEF_1:
    if i['gene_name'] in union_gene['gene_name']:
        gene_id.append((i['geneID'] , i['gene_name']))

gene_id = set(gene_id)
gene_id = list(gene_id)       


#for i in MEF_2:
#    if i['gene_name'] in union_gene['gene_name']:
#        gene_id.append((i['geneID'] , i['gene_name']))
#        
#gene_id = set(gene_id)
#gene_id = list(gene_id)       
#
#
#for i in IPS_P3_1:
#    if i['gene_name'] in union_gene['gene_name']:
#        gene_id.append((i['geneID'] , i['gene_name']))
#        
#gene_id = set(gene_id)
#gene_id = list(gene_id)    
#
#for i in IPS_P3_2:
#    if i['gene_name'] in union_gene['gene_name']:
#        gene_id.append((i['geneID'] , i['gene_name']))
#        
#gene_id = set(gene_id)
#gene_id = list(gene_id)    
#
#
#for i in E14_1:
#    if i['gene_name'] in union_gene['gene_name']:
#        gene_id.append((i['geneID'] , i['gene_name']))
#        
#gene_id = set(gene_id)
#gene_id = list(gene_id)   
#
#
#for i in E14_2:
#    if i['gene_name'] in union_gene['gene_name']:
#        gene_id.append((i['geneID'] , i['gene_name']))
#        
#gene_id = set(gene_id)
#gene_id = list(gene_id)   

gene_id = np.array(gene_id , dtype = data_type)

MEF_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\MEF_1_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)
MEF_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\MEF_2_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)
IPS_P3_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\IPS_P3_1_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)
IPS_P3_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\IPS_P3_2_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)
E14_1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\E14_1_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)
E14_2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\RNA\\E14_2_normalized_signal.txt' , dtype = data_type_1 , skiprows = 5)



IPSC_gene_expression = []
for i in set(gene_id['gene_name']):
    gene = union_gene[union_gene['gene_name'] == i][0]
    gene_ID = gene_id[gene_id['gene_name'] == i][0]['geneID']
    MEF_FPKM_1 = MEF_1[MEF_1['geneID'] == gene_ID]['FPKM'].mean()
    MEF_FPKM_2 = MEF_2[MEF_2['geneID'] == gene_ID]['FPKM'].mean()
    IPS_P3_FPKM_1 = IPS_P3_1[IPS_P3_1['geneID'] == gene_ID]['FPKM'].mean()
    IPS_P3_FPKM_2 = IPS_P3_2[IPS_P3_2['geneID'] == gene_ID]['FPKM'].mean()
    E14_FPKM_1 = E14_1[E14_1['geneID'] == gene_ID]['FPKM'].mean()
    E14_FPKM_2 = E14_2[E14_2['geneID'] == gene_ID]['FPKM'].mean()
    
    IPSC_gene_expression.append((gene['gene_name']  , gene['chr'] , gene['strand'] , gene['start'] , gene['end'] , gene['CCS'] , gene['NT5'] , gene['NT6'] , gene['fESC'] , (MEF_FPKM_1 + MEF_FPKM_2) / 2 , (IPS_P3_FPKM_1 + IPS_P3_FPKM_2) / 2 , (E14_FPKM_1 + E14_FPKM_2) / 2))

IPSC_gene_expression = np.array(IPSC_gene_expression , dtype = union_type_1)


for i in range(len(np.isnan(IPSC_gene_expression['MEF']))):
    if np.isnan(IPSC_gene_expression['MEF'])[i] == True:
        print IPSC_gene_expression[i]
        IPSC_gene_expression[i]['MEF'] = 0
        IPSC_gene_expression[i]['IPS_P3'] = 0
        IPSC_gene_expression[i]['E14'] = 0
        print IPSC_gene_expression[i]
        
        
Write2fils_nochr('H:\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt' , IPSC_gene_expression)

