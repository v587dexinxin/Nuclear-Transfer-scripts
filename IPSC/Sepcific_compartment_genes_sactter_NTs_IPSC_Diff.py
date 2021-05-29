# -*- coding: utf-8 -*-
"""
Created on Tue May 11 22:02:02 2021

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
# import cPickle
import pandas as pd

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')



    
union_gene = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\reads_count\\Filtered_NTs_IPSC_diff_genes_p_0.05.csv' , index_col = 1)
union_gene = pd.concat([union_gene['Gene_ID'] , pd.DataFrame(union_gene.index , index=union_gene.index) , union_gene['Chr'] , union_gene['Strand'] , 
                        union_gene['Start'], union_gene['End'], union_gene['baseMean'], union_gene['log2FoldChange'], 
                        union_gene['lfcSE'], union_gene['stat'], union_gene['pvalue'], union_gene['padj'], 
                        union_gene[['NT5_1' , 'NT5_2' , 'NT5_3' , 'NT5_4' , 'NT6_1' , 'NT6_2' , 'NT6_3']].mean(1) , 
                        union_gene[['IPSC_1' , 'IPSC_2']].mean(1)] , axis = 1)

union_gene.columns = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'baseMean' , 'log2FoldChange' , 'lfcSE' ,  
                      'stat' , 'pvalue' , 'padj' , 'NTs_FPKM' , 'IPSC_FPKM']



cl1 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_A_B_Method_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl2 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_A_B_Donor_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl3 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_A_B_Method_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl4 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_A_B_Donor_specific_genes.txt' , dtype = 'U64' , usecols = (0))

cl5 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_B_A_Method_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl6 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\NT_B_A_Donor_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl7 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_B_A_Method_specific_genes.txt' , dtype = 'U64' , usecols = (0))
cl8 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_Compartment_related_genes\\IPSC_B_A_Donor_specific_genes.txt' , dtype = 'U64' , usecols = (0))



diff = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\reads_count\\Filtered_NTs_IPSC_diff_genes_q_0.05.csv' , index_col=1)

cl1_gene = [] ; cl2_gene = [] ; cl3_gene = [] ; cl4_gene = []
cl5_gene = [] ; cl6_gene = [] ; cl7_gene = [] ; cl8_gene = []


for i in cl1:
    if i in diff.index:
        cl1_gene.append(i)
        
        
for i in cl2:
    if i in diff.index:
        cl2_gene.append(i)
        
        
for i in cl3:
    if i in diff.index:
        cl3_gene.append(i)
        
for i in cl4:
    if i in diff.index:
        cl4_gene.append(i)

for i in cl5:
    if i in diff.index:
        cl5_gene.append(i)
        
for i in cl6:
    if i in diff.index:
        cl6_gene.append(i)

for i in cl7:
    if i in diff.index:
        cl7_gene.append(i)

for i in cl8:
    if i in diff.index:
        cl8_gene.append(i)        
        
out = pd.ExcelWriter('Specific_NTs_IPSC_diff_genes.xlsx')
union_gene.loc[cl1_gene,:].to_excel(out, sheet_name='NT_A_B_Method_genes' , index = 0)
union_gene.loc[cl2_gene,:].to_excel(out, sheet_name='NT_A_B_Donor_genes' , index = 0)
union_gene.loc[cl3_gene,:].to_excel(out, sheet_name='IPSC_A_B_Method_genes' , index = 0)
union_gene.loc[cl4_gene,:].to_excel(out, sheet_name='IPSC_A_B_Donor_genes' , index = 0)
union_gene.loc[cl5_gene,:].to_excel(out, sheet_name='NT_B_A_Method_genes' , index = 0)
union_gene.loc[cl6_gene,:].to_excel(out, sheet_name='NT_B_A_Donor_genes' , index = 0)
union_gene.loc[cl7_gene,:].to_excel(out, sheet_name='IPSC_B_A_Method_genes' , index = 0)
union_gene.loc[cl8_gene,:].to_excel(out, sheet_name='IPSC_B_A_Donor_genes' , index = 0)

out.save()




pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\Specific_compartment_related_genes\\Specific_A_B_genes.pdf')
size = (10, 10)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
fig = plt.figure(figsize = size)
ax = fig.add_axes([Left , HB , width, HH])
# ax.scatter(np.log2(np.array(union_gene.loc[cl2_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl2_gene].NTs_FPKM) + 1), alpha = 1 , c = 'blue')
ax.scatter(np.log2(np.array(union_gene.loc[cl1_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl1_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'deeppink')
# ax.scatter(np.log2(np.array(union_gene.loc[cl6_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl6_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'lightskyblue')
ax.scatter(np.log2(np.array(union_gene.loc[cl5_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl5_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'sandybrown')

ax.set_xlabel('IPSC_log2(FPKM+1)' , size = 20)
ax.set_ylabel('NTs_log2(FPKM+1)' , size = 20)

ax.set_xlim(-0.5 , 12)
ax.set_ylim(-0.5 , 12)

ax.plot([8 , 9] , [9.5 , 9.5] , c = 'sandybrown' , lw = 2.5 )
# ax.plot([8 , 9] , [10 , 10] , c = 'lightskyblue' , lw = 2.5 )
ax.plot([8 , 9] , [10.5 , 10.5] , c = 'deeppink' , lw = 2.5 )
# ax.plot([8 , 9] , [11 , 11] , c = 'blue' , lw = 2.5 )

ax.text(9.3 , 9.4 , 'NT_B_A_Method' , size = 10 )
# ax.text(9.3 , 9.9 , 'NT_B_A_Donor' , size = 10 )
ax.text(9.3 , 10.4 , 'NT_A_B_Method' , size = 10 )
# ax.text(9.3 , 10.9, 'NT_A_B_Donor' , size = 10 )

ax.set_title('NT process' , size = 20)

pp.savefig(fig)
pp.close()       


pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\Specific_compartment_related_genes\\Specific_B_A_genes.pdf')
size = (10, 10)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
fig = plt.figure(figsize = size)
ax = fig.add_axes([Left , HB , width, HH])
# ax.scatter(np.log2(np.array(union_gene.loc[cl4_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl4_gene].NTs_FPKM) + 1), alpha = 1 , c = 'blue')
ax.scatter(np.log2(np.array(union_gene.loc[cl3_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl3_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'deeppink')
# ax.scatter(np.log2(np.array(union_gene.loc[cl8_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl8_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'lightskyblue')
ax.scatter(np.log2(np.array(union_gene.loc[cl7_gene].IPSC_FPKM) + 1) , np.log2(np.array(union_gene.loc[cl7_gene].NTs_FPKM) + 1) , alpha = 1 , c = 'sandybrown')

ax.set_xlabel('IPSC_log2(FPKM+1)' , size = 20)
ax.set_ylabel('NTs_log2(FPKM+1)' , size = 20)

ax.set_xlim(-0.5 , 12)
ax.set_ylim(-0.5 , 12)

ax.plot([8 , 9] , [9.5 , 9.5] , c = 'sandybrown' , lw = 2.5 )
# ax.plot([8 , 9] , [10 , 10] , c = 'lightskyblue' , lw = 2.5 )
ax.plot([8 , 9] , [10.5 , 10.5] , c = 'deeppink' , lw = 2.5 )
# ax.plot([8 , 9] , [11 , 11] , c = 'blue' , lw = 2.5 )

ax.text(9.3 , 9.4 , 'IPSC_B_A_Method' , size = 10 )
# ax.text(9.3 , 9.9 , 'IPSC_B_A_Donor' , size = 10 )
ax.text(9.3 , 10.4 , 'IPSC_A_B_Method' , size = 10 )
# ax.text(9.3 , 10.9, 'IPSC_A_B_Donor' , size = 10 )

ax.set_title('IPSC process' , size = 20)

pp.savefig(fig)
pp.close()       




