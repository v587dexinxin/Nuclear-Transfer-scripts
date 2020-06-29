# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:36:26 2020

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

my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})

locus_type = np.dtype({'names':['chr' , 'locus'] , 
                      'formats':['S8' , np.int]})

def Get_specificAB(compartment1 , compartment2):
    specific_CCS_A = []
    specific_CCS_B = []
    specific_MEF_A = []
    specific_MEF_B = []
    for g in chro:
        tmp1 = compartment1[compartment1['chr'] == g]
        tmp2 = compartment2[compartment2['chr'] == g]
        for i in range(len(tmp1)):
            site = int(i * res)
            if tmp1[i]['pc'] > 0 and tmp2[i]['pc'] < 0:
                specific_CCS_A.append((g , site))
                specific_MEF_B.append((g , site))
            elif tmp1[i]['pc'] < 0 and tmp2[i]['pc'] > 0:
                specific_CCS_B.append((g , site))
                specific_MEF_A.append((g , site))
            else:
                pass
    specific_CCS_A = np.array(specific_CCS_A , dtype = locus_type)
    specific_CCS_B = np.array(specific_CCS_B , dtype = locus_type)
    specific_MEF_A = np.array(specific_MEF_A , dtype = locus_type)
    specific_MEF_B = np.array(specific_MEF_B , dtype = locus_type)
    return specific_CCS_A , specific_CCS_B , specific_MEF_A , specific_MEF_B
    
    
def Common2(comp1 , comp2):
    common = []
    for g in chro:
        tmp1 = comp1[comp1['chr'] == g]
        tmp2 = comp2[comp2['chr'] == g]
        for i in tmp1:
            if i['locus'] in tmp2['locus']:
                common.append(i)
    common = np.array(common , dtype = locus_type)
    return common
    
def Common3(comp1 , comp2 , comp3):
    common = []
    for g in chro:
        tmp1 = comp1[comp1['chr'] == g]
        tmp2 = comp2[comp2['chr'] == g]
        tmp3 = comp3[comp3['chr'] == g]
        for i in tmp1:
            if (i['locus'] in tmp2['locus']) and (i['locus'] in tmp3['locus']):
                common.append(i)
    common = np.array(common , dtype = locus_type)
    return common

       
def Common4(comp1 , comp2 , comp3 , comp4):
    common = []
    for g in chro:
        tmp1 = comp1[comp1['chr'] == g]
        tmp2 = comp2[comp2['chr'] == g]
        tmp3 = comp3[comp3['chr'] == g]
        tmp4 = comp4[comp4['chr'] == g]
        for i in tmp1:
            if (i['locus'] in tmp2['locus']) and (i['locus'] in tmp3['locus']) and (i['locus'] in tmp4['locus']):
                common.append(i)
    common = np.array(common , dtype = locus_type)
    return common 


def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})
    site_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union_gene_1 = []
    for i in union_gene:
        union_gene_1.append((i['gene_name'] , i['chr'] , (i['start'] + i['end']) / 2 , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
    union_gene = np.array(union_gene_1 , dtype = site_type)
    return union_gene
    
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
    
    
def Get_Venn4_data(data_set , classify):
    data = {classify[0] : data_set[0] , classify[1] : data_set[1] , 
            classify[2] : data_set[2] , classify[3] : data_set[3]}
                
            
    common_12 = Common2(data[classify[0]] , data[classify[1]])
    common_13 = Common2(data[classify[0]] , data[classify[2]])
    common_14 = Common2(data[classify[0]] , data[classify[3]])
    common_23 = Common2(data[classify[1]] , data[classify[2]])
    common_24 = Common2(data[classify[1]] , data[classify[3]])
    common_34 = Common2(data[classify[2]] , data[classify[3]])
    common_123 = Common3(data[classify[0]] , data[classify[1]] , data[classify[2]])
    common_124 = Common3(data[classify[0]] , data[classify[1]] , data[classify[3]])
    common_134 = Common3(data[classify[0]] , data[classify[2]] , data[classify[3]])
    common_234 = Common3(data[classify[1]] , data[classify[2]] , data[classify[3]])
    common_1234 = Common4(data[classify[0]] , data[classify[1]] , data[classify[2]] , data[classify[3]])
    
    
    
    print 'area1=' + str(len(data[classify[0]]))
    print 'area2=' + str(len(data[classify[1]]))
    print 'area3=' + str(len(data[classify[2]]))
    print 'area4=' + str(len(data[classify[3]]))
    print 'n12=' + str(len(common_12)) 
    print 'n13=' + str(len(common_13))
    print 'n14=' + str(len(common_14))
    print 'n23=' + str(len(common_23))
    print 'n24=' + str(len(common_24))
    print 'n34=' + str(len(common_34))
    print 'n123=' + str(len(common_123))
    print 'n124=' + str(len(common_124))
    print 'n134=' + str(len(common_134))
    print 'n234=' + str(len(common_234))
    print 'n1234=' + str(len(common_1234))
    return data 
    
def Write_Venn4_partial_genes(data , classify , dynamic):
    data1 = []
    data2 = []
    data3 = []
    data4 = []
    data5 = []
    data6 = []
    a = dynamic.split('_')
    
    
    for g in chro:
        tmp1 = data[classify[0]][data[classify[0]]['chr'] == g]
        tmp2 = data[classify[2]][data[classify[2]]['chr'] == g]
        for i in tmp1:
            if i['locus'] not in tmp2['locus']:
                data1.append(i)
            else:
                data2.append(i)
        for i in tmp2:
            if i['locus'] not in tmp1['locus']:
                data3.append(i)
    
    data1 = np.array(data1 , dtype = locus_type)
    data2 = np.array(data2 , dtype = locus_type)
    data3 = np.array(data3 , dtype = locus_type)
        
    for g in chro:
        tmp1 = data[classify[1]][data[classify[1]]['chr'] == g]
        tmp2 = data[classify[3]][data[classify[3]]['chr'] == g]
        for i in tmp1:
            if i['locus'] not in tmp2['locus']:
                data6.append(i)
            else:
                data5.append(i)
        for i in tmp2:
            if i['locus'] not in tmp1['locus']:
                data4.append(i)
                
    data4 = np.array(data4 , dtype = locus_type)
    data5 = np.array(data5 , dtype = locus_type)
    data6 = np.array(data6 , dtype = locus_type)       
    
        
    data = {'CCS_' + dynamic[0] : data1 , 'CCS_' + dynamic[0] + '_NT_'+ dynamic[0] + '_' + dynamic[1] : data2 , 'NT_' + dynamic[0] + '_' + dynamic[1] : data3 , 
            'IPSC_' + dynamic[0] + '_' + dynamic[1] : data4 , 'IPSC_' + dynamic[0] + '_' + dynamic[1] + '_MEF_' + dynamic[0] : data5 , 'MEF_' + dynamic[0] : data6}
            
            
    for k in data:
        if len(a) > 1:
            out = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_resist_' + dynamic[0] + '_' + dynamic[1] + '\\Resist_' + k + '_genes.txt', 'w')
        else:
            out = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_' + dynamic[0] + '_' + dynamic[1] + '\\' + k + '_genes.txt', 'w')
        gene = get_compartment_genes(data[k] , union_gene)
        print k , len(data[k]) , len(gene)
        for i in gene:
            out.writelines(i['gene_name'] + '\n')
        out.close()
        


    
chro = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']                   
res = 200000
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']

MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']

specific_NT_A_B = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Fig3E_dynamic_genes_NT_vs_IPSC\\NT_noIPSC_A_B.txt' , dtype = locus_type)
specific_NT_B_A = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Fig3E_dynamic_genes_NT_vs_IPSC\\NT_noIPSC_B_A.txt' , dtype = locus_type)
specific_IPSC_A_B = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Fig3E_dynamic_genes_NT_vs_IPSC\\noNT_IPSC_A_B.txt' , dtype = locus_type)
specific_IPSC_B_A = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Fig3E_dynamic_genes_NT_vs_IPSC\\noNT_IPSC_B_A.txt' , dtype = locus_type)


union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

specific_CCS_A , specific_CCS_B , specific_MEF_A , specific_MEF_B = Get_specificAB(CCS , MEF)

##A_B------------------------------


classify = ['specific_CCS_A' , 'specific_MEF_A' , 'specific_NT_A_B' , 'specific_IPSC_A_B']
data_set = [specific_CCS_A , specific_MEF_A , specific_NT_A_B , specific_IPSC_A_B]
    
data = Get_Venn4_data(data_set , classify)
Write_Venn4_partial_genes(data , classify , 'AB')


##B_A------------------------------


classify = ['specific_CCS_B' , 'specific_MEF_B' , 'specific_NT_B_A' , 'specific_IPSC_B_A']
data_set = [specific_CCS_B , specific_MEF_B , specific_NT_B_A , specific_IPSC_B_A]
    
data = Get_Venn4_data(data_set , classify)
Write_Venn4_partial_genes(data , classify , 'BA')





##---------------------------Resist-----------------------------------------
   
specific_NT_resist_A_B = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Respective_resist_compartment\\NT_noIPSC_resist_A_B.txt' , dtype = locus_type)
specific_NT_resist_B_A = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Respective_resist_compartment\\NT_noIPSC_resist_B_A.txt' , dtype = locus_type)
specific_IPSC_resist_A_B = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Respective_resist_compartment\\noNT_IPSC_resist_A_B.txt' , dtype = locus_type)
specific_IPSC_resist_B_A = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Respective_resist_compartment\\noNT_IPSC_resist_B_A.txt' , dtype = locus_type)


specific_CCS_A , specific_CCS_B , specific_MEF_A , specific_MEF_B = Get_specificAB(CCS , MEF)


##A_B------------------------------
data_set = [specific_CCS_A , specific_MEF_A , specific_NT_resist_A_B , specific_IPSC_resist_A_B]
classify = ['specific_CCS_A' , 'specific_MEF_A' , 'specific_NT_resist_A_B' , 'specific_IPSC_resist_A_B']

    
data = Get_Venn4_data(data_set , classify)
Write_Venn4_partial_genes(data , classify , 'AB_resist')



##B_A------------------------------
        
data_set = [specific_CCS_B , specific_MEF_B , specific_NT_resist_B_A , specific_IPSC_resist_B_A]
classify = ['specific_CCS_B' , 'specific_MEF_B' , 'specific_NT_resist_B_A' , 'specific_IPSC_resist_B_A']

    
data = Get_Venn4_data(data_set , classify)     
Write_Venn4_partial_genes(data , classify , 'BA_resist')
   

        
