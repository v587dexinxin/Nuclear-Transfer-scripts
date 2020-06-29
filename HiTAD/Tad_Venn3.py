# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:15:39 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})


union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
         
                 
                 

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()
    
def Load_Tads(TadFil):
    Tads = np.loadtxt(TadFil , dtype = tad_type)
    return Tads

def union_Tads(data_list):
    Tads_all = []
    for c in data_list:
        for i in c:
            Tads_all.append(i)
    Tads_all.sort(key = lambda x:(x[0],x[1],x[2]))
    Tads_all = np.array(Tads_all ,dtype = tad_type)
    union_tads = []
    for g in chrom:
        Tads = Tads_all[Tads_all['chr'] == g]
        Tads_number = len(Tads)
        for i in range(Tads_number-1 , -1 , -1):
            start = Tads[i]['start']
            end = Tads[i]['end']
            prev_ss = Tads[i-1]['start'] - 40000
            prev_se = Tads[i-1]['start'] + 40000
            prev_es = Tads[i-1]['end'] - 40000
            prev_ee = Tads[i-1]['end'] + 40000
            if (start >= prev_ss) and (start <= prev_se) and (end >= prev_es) and (end <= prev_ee):
                Tads = list(Tads)
                Tads.remove(Tads[i])
                Tads = np.array(Tads,dtype = tad_type)
            else:
                pass
        union_tads.extend(list(Tads))
            
    union_tads = np.array(union_tads , dtype = tad_type) 
    return union_tads
    
    
    
def Common_Tads(data_list , union_tads):
    commom = []
    for i in union_tads:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        n = 0
        for j in data_list:
            tmp = j[j['chr'] == chro]
            mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
            overlap = tmp[mask]
            if overlap.size != 0:
                n += 1
            else:
                pass
        if n == len(data_list):
            commom.append(i)
    commom = np.array(commom , dtype = tad_type)
    return commom

def Common2_Tads(tads1 , tads2):
    common = []
    remain1 = []
    remain2 = []
    n1 = 0
    n2 = 0
    for i in tads1:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tads2[tads2['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size != 0:
            common.append(i)
            n1 += 1
        else:
            remain1.append(i)
    for i in tads2:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tads1[tads1['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size != 0:
            n2 += 1
        else:
            remain2.append(i)
    common = np.array(common , dtype = tad_type)
    remain1 = np.array(remain1 , dtype = tad_type)
    remain2 = np.array(remain2 , dtype = tad_type)
    
    return common , remain1 , remain2
    
def diff_Tads(tads1 , tads2):
    diff1 = []
    common = []
    for i in tads1:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tads2[tads2['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size == 0:
            diff1.append(i)
        else:
            common.append(i)
    diff1 = np.array(diff1 , dtype = tad_type)
    common = np.array(common , dtype = tad_type)
    return diff1 , common
    
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
        
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
        for i in gene:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()                        
            
TadFolder = 'H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19']
cells = ['NT5','NT6','fESC']

CCS = Load_Tads(os.path.join(TadFolder , 'CCS_40K_allreps_union_small_domain.txt'))
CCS = CCS[CCS['chr'] != 'X']
NT5 = Load_Tads(os.path.join(TadFolder , 'NT5_40K_allreps_union_small_domain.txt'))
NT5 = NT5[NT5['chr'] != 'X']
NT6 = Load_Tads(os.path.join(TadFolder , 'NT6_40K_allreps_union_small_domain.txt'))
NT6 = NT6[NT6['chr'] != 'X']
fESC = Load_Tads(os.path.join(TadFolder , 'fESC_40K_allreps_union_small_domain.txt'))
fESC = fESC[fESC['chr'] != 'X']



common_NT , n1 , n2= Common2_Tads(NT5 , NT6)
union_CCS_NT_fESC = union_Tads([CCS , common_NT , fESC])

CCS_NT , CCS_noNT , noCCS_NT = Common2_Tads(CCS , common_NT)
CCS_fESC , CCS_nofESC , noCCS_fESC = Common2_Tads(CCS , fESC)
NT_fESC , NT_nofESC , noNT_fESC = Common2_Tads(common_NT , fESC)
CCS_NT_fESC = Common_Tads([CCS , common_NT , fESC] , union_CCS_NT_fESC)


CCS_noNT_nofESC , CCS_noNT_fESC= diff_Tads(CCS_noNT , fESC)
CCS_nofESC_noNT , CCS_nofESC_NT= diff_Tads(CCS_nofESC , common_NT)
NT_nofESC_noCCS , NT_nofESC_CCS = diff_Tads(NT_nofESC , CCS)

noCCS_NT_nofESC , noCCS_NT_fESC = diff_Tads(noCCS_NT , fESC)
noCCS_fESC_noNT , noCCS_fESC_NT = diff_Tads(noCCS_fESC , common_NT)
noNT_fESC_noCCS , noNT_fESC_CCS = diff_Tads(noNT_fESC , CCS)


Tads = {'CCS_noNT_nofESC' : CCS_noNT_nofESC , 'CCS_noNT_fESC' : CCS_noNT_fESC , 'CCS_NT_nofESC' : CCS_nofESC_NT,
        'noCCS_NT_fESC' : noCCS_fESC_NT , 'noCCS_NT_nofESC' : NT_nofESC_noCCS , 'noCCS_noNT_fESC' : noCCS_fESC_noNT }
        
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)



genes = {}

for cl in Tads:
    genes[cl] = []
    for i in Tads[cl]:
        chro = i['chr']
        start = i['start']
        end = i['end']
        tmp_gene = union[union['chr'] == 'chr' + chro]        
        mask = (tmp_gene['gene_site'] >= start) & (tmp_gene['gene_site'] <= end)
        overlap = tmp_gene[mask]
        if overlap.size != 0 :
            for j in overlap:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
                
for k,v in genes.items():
    Write2genes('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\TAD_Venn3_genes\\' + k + '_genes.txt' , v)                


Abc = len(CCS) - len(CCS_NT) - len(CCS_fESC) + len(CCS_NT_fESC)
aBc = len(common_NT) - len(CCS_NT) - len(NT_fESC) + len(CCS_NT_fESC)
ABc = len(CCS_NT) - len(CCS_NT_fESC)
abC = len(fESC) - len(CCS_fESC) - len(NT_fESC) + len(CCS_NT_fESC)
AbC = len(CCS_fESC) - len(CCS_NT_fESC)
aBC = len(NT_fESC) - len(CCS_NT_fESC)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC


fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('CCS', 'NT' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig5_figs\\CCS_NT(common)fESC_Venn3.pdf')


Ab = len(CCS) - len(CCS_fESC)
aB = len(fESC) - len(CCS_fESC)
AB = len(CCS_fESC)

fig = plt.figure(figsize = (10, 10))
venn2(subsets = (Ab , aB , AB), set_labels = ('CCS', 'fESC'))

run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig5_figs\\CCS_fESC_Venn2.pdf')




union_NT = union_Tads([NT5 , NT6])
union_CCS_NT_fESC = union_Tads([CCS , union_NT , fESC])


CCS_NT , n1 , n2 = Common2_Tads(CCS , union_NT)
CCS_fESC , n1 , n2 = Common2_Tads(CCS , fESC)
NT_fESC , n1 , n2 = Common2_Tads(union_NT , fESC)
CCS_NT_fESC = Common_Tads([CCS , union_NT , fESC] , union_CCS_NT_fESC)

Abc = len(CCS) - len(CCS_NT) - len(CCS_fESC) + len(CCS_NT_fESC)
aBc = len(union_NT) - len(CCS_NT) - len(NT_fESC) + len(CCS_NT_fESC)
ABc = len(CCS_NT) - len(CCS_NT_fESC)
abC = len(fESC) - len(CCS_fESC) - len(NT_fESC) + len(CCS_NT_fESC)
AbC = len(CCS_fESC) - len(CCS_NT_fESC)
aBC = len(NT_fESC) - len(CCS_NT_fESC)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC

fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('CCS', 'NT' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\CCS_NT(union)fESC_Venn3.pdf')



  