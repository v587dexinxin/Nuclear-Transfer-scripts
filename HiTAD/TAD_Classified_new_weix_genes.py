# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 02:09:50 2021

@author: xxli
"""

import numpy as np 
from itertools import islice
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_union_gene_sites(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if (i['CCS_FPKM'] > 0.1) or (i['NT5_FPKM'] > 0.1) or (i['NT6_FPKM'] > 0.1) or (i['F35_FPKM'] > 0.1) or (i['F40_FPKM'] > 0.1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS_FPKM'] , i['NT2_FPKM'] , i['NT3_FPKM'] , i['NT4_FPKM'] , i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F37_FPKM'] , i['F40_FPKM'] , i['F41_FPKM']))
    union = np.array(union , dtype = union_type_1)
    return union

def get_raw_genes_new(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT2' , 'NT3' , 'NT4' , 'NT5' , 'NT6' , 'F35' , 'F37' , 'F40' , 'F41'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene


def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int , np.int]})
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

def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'log2FoldChange']) + '\n')
        for i in gene:
            out.writelines('\t'.join([i['gene_name'] , i['chr'] , str(i['gene_site']) , str(i['CCS_FPKM']) , str(i['NT5_FPKM']) , str(i['NT6_FPKM']) , str(i['F35_FPKM']) , str(i['F40_FPKM']) , str(i['log2FoldChange'])])+'\n')
    out.close()  

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            out.writelines(i +'\n')
    out.close()
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
   
    
def OutputTXT(Info , outfil , head):
    """
    """
    print ("Output ...")
    
    # head = ['chr','position','PC-M','PC-P','diff','P_Value']
    with open(outfil, 'w') as out:
        out.writelines('\t'.join(head)+'\n')
        
        for line in Info:
            line = map(str, line)
            out.writelines('\t'.join(line)+'\n')
    out.close()
    

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['U8' , np.int , np.int]})


Repro = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\Repro_TADs.txt' , dtype = tad_type )
Resis = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\Resis_TADs.txt' , dtype = tad_type )
Over = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\Over_TADs.txt' , dtype = tad_type )
Others = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\Others_TADs.txt' , dtype = tad_type )

A_B_Repro = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\A_B_Repro_TADs.txt' , dtype = tad_type )
A_B_Resis = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\A_B_Resis_TADs.txt' , dtype = tad_type )
A_B_Over = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\A_B_Over_TADs.txt' , dtype = tad_type )

B_A_Repro = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\B_A_Repro_TADs.txt' , dtype = tad_type )
B_A_Resis = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\B_A_Resis_TADs.txt' , dtype = tad_type )
B_A_Over = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\B_A_Over_TADs.txt' , dtype = tad_type )



tads = {'Repro':Repro , 'Resis':Resis , 'Over':Over , 'Others':Others , 'A_B_Repro':A_B_Repro , 'A_B_Resis':A_B_Resis,
        'A_B_Over':A_B_Over , 'B_A_Repro':B_A_Repro , 'B_A_Resis':B_A_Resis , 'B_A_Over':B_A_Over}



gtf = Load_gtf('E:\\Data\\literature_data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene_sites('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')
union_gene_raw = get_raw_genes_new('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')

genes = {}
gene_id = {}


for cl in tads:
    genes[cl] = []
    gene_id[cl] = []
    for i in tads[cl]:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000 
        tmp_gene = union_gene[union_gene['chr'] == 'chr' + chro]    
        mask1 = (tmp_gene['gene_site'] >= start_s) & (tmp_gene['gene_site'] <= start_e)
        mask2 = (tmp_gene['gene_site'] >= end_s) & (tmp_gene['gene_site'] <= end_e)
        overlap1 = tmp_gene[mask1]
        overlap2 = tmp_gene[mask2]
        if overlap1.size != 0 :
            for j in overlap1:
                tmp_gtf = gtf[gtf['gene_name'] == j['gene_name']][0]
                genes[cl].append(j)
                gene_id[cl].append(tmp_gtf['gene_id'])
        if overlap2.size != 0 :
            for j in overlap2:
                tmp_gtf = gtf[gtf['gene_name'] == j['gene_name']][0]
                genes[cl].append(j)
                gene_id[cl].append(tmp_gtf['gene_id'])
    genes[cl] = np.array(genes[cl] , dtype = union_gene.dtype)
     
    
    
head = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCs_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F40_FPKM']
   
for k,v in genes.items():
    tmp_genes = [] ; tmp_gene_id = []      
    print (k , len(v))
    gene_name = set(v['gene_name'])
    Gene = []
    for i in gene_name:
        gene = union_gene[union_gene['gene_name'] == i][0]
        gene1 = union_gene_raw[union_gene_raw['gene_name'] == i][0]
        g_id = gtf[gtf['gene_name'] == i][0]['gene_id']
        tmp_genes.append(gene)
        tmp_gene_id.append(g_id)
        Gene.append([g_id , gene1['gene_name'] , gene1['chr'] , gene1['strand'] , gene1['start'] , gene1['end'] , gene['CCS_FPKM'] , gene['NT5_FPKM'] , gene['NT6_FPKM'] , gene['F35_FPKM'] , gene['F40_FPKM']])
    tmp_genes = np.array(tmp_genes)
    print (k , len(tmp_genes))
    # Write2genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\TAD_related_genes\\' + k + '_genes.txt' , tmp_genes)                
    # Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\TAD_related_genes\\' + k + '_genes_ID.txt' , tmp_gene_id)       
    OutputTXT(Gene , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\TAD_related_genes\\excel\\' + k + '_genes.txt' , head)


##--------------------------DEG-TAD_percent----------------------------




tad = {'Repro':Repro , 'Resis':Resis , 'Over':Over , 'Others':Others}

Diff = pd.read_csv(open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\diff_expression_q_0.05\\RNA_fESC_vs_NTs_diff_expression_q_0.05.csv' , 'r') , index_col = 0)
gene_names = []
for i in Diff.index:
    g_name = gtf[gtf['gene_id'] == i][0]['gene_name']
    if g_name == 'U2af1l4':
        print (i)
    gene_names.append(g_name)
    
diff_gene = set(gene_names)


union_diff = []
for i in diff_gene:
    overlap = union_gene[union_gene['gene_name'] == i][0]
    union_diff.append(overlap)
union_diff = np.array(union_diff ,  dtype=union_gene.dtype)
    
  
percent = {}
genes_diff = {}
for cl in tad:
    n = 0
    genes_diff[cl] = []
    for i in tad[cl]:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = union_diff[union_diff['chr'] == 'chr' + chro]
        mask1 = (tmp['gene_site'] >= start_s) & (tmp['gene_site'] <= start_e)
        mask2 = (tmp['gene_site'] >= end_s) & (tmp['gene_site'] <= end_e)
        
        overlap1 = tmp[mask1]
        overlap2 = tmp[mask2]
        if overlap1.size != 0:
            for j in overlap1:
                if j['gene_name'] == 'U2af1l4':
                    print (j)
                genes_diff[cl].append(j)
        if overlap2.size != 0:
            for j in overlap2:
                if j['gene_name'] == 'U2af1l4':
                    print (j)
                genes_diff[cl].append(j)
        if (overlap1.size != 0) or (overlap2.size != 0):
            n += 1
    genes_diff[cl] = np.array(genes_diff[cl] , dtype = union_gene.dtype)
    percent[cl] = n / len(tad[cl])
    print (cl , n / len(tad[cl]))



  
diff = {}  
# percent = {}
for cl in genes_diff:    
    diff[cl] = []
    for i in set(genes_diff[cl]['gene_name']):
        if i in diff_gene:
            tmp = genes_diff[cl][genes_diff[cl]['gene_name'] == i][0]
            diff[cl] .append(tmp)
    diff[cl] = np.array(diff[cl] , dtype = union_gene.dtype)
    # percent_b = len(diff[cl]) / len(genes[cl])
    # percent[cl] = percent_b
    # print (cl , percent_b)



fig = plt.figure(figsize = (10, 10))
ax = fig.add_axes([0.2 , 0.2 , 0.6 , 0.6])
ax.bar([1 , 2, 3 , 4] , [percent['Repro'] , percent['Resis'] , percent['Over'] , percent['Others']] , color=['firebrick' , 'green' , 'skyblue' , 'plum'])
ax.text(0.7 , percent['Repro'] + 0.05 ,  np.round(percent['Repro'] , 4))
ax.text(1.7 , percent['Resis'] + 0.05 ,  np.round(percent['Resis'] , 4))
ax.text(2.7 , percent['Over'] + 0.05 ,  np.round(percent['Over'] , 4))
ax.text(3.7 , percent['Others'] + 0.05 ,  np.round(percent['Others'] , 4))

ax.set_xticks([1 , 2 , 3 , 4])
ax.set_xticklabels(['Repro' , 'Resis' , 'Over' , 'Static'])
ax.set_xlim((0.5 , 4.5))
ax.set_ylim((0 , 0.4))
ax.set_ylabel('DEG percent')

run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\plot\\Tad_DEG_percent.pdf')


##----------------------------Write_to_DEG_Files----------------------------------------------

union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM' , 'log2FoldChange'],
             'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
     
    
for k,v in diff.items():
    tmp_genes = [] ; tmp_gene_id = []      
    print (k , len(v))
    gene_name = set(v['gene_name'])
    for i in gene_name:
        
        g_ids = gtf[gtf['gene_name'] == i]['gene_id']
        for j in g_ids:
            if j in Diff.index:
                g_id = j
                        
        gene = tuple(list(union_gene[union_gene['gene_name'] == i][0]) + [Diff.loc[g_id]['log2FoldChange']])
        tmp_genes.append(gene)
        tmp_gene_id.append(g_id)
    tmp_genes = np.array(tmp_genes , dtype = union_type_1)
    print (k , len(tmp_genes))
    Write2genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\TAD_related_genes\\TAD_related_DEGs\\' + k + '_DEGs.txt' , tmp_genes)                
    Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\TAD_classify_weix\\TAD_related_genes\\TAD_related_DEGs\\' + k + '_DEGs_ID.txt' , tmp_gene_id)       

















