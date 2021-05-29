# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 14:43:59 2020

@author: xxli
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Our Own Color Map
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')


d_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35_R1' , 'F35_R2' , 'F35_R3' , 'F37_R1' , 'F37_R2' , 'F37_R3' , 'F41_R1' , 'F41_R2' , 'F41_R3' , 'NT2_R1' , 'NT2_R2' , 'NT2_R3' , 'NT3_R1' , 'NT3_R2' , 'NT3_R3' , 'NT4_R1' , 'NT4_R2' , 'NT4_R3'],
                   'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


d_type_1 = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35' , 'F37' , 'F41' , 'NT2' , 'NT3' , 'NT4'],
                     'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float]})

                

def get_raw_genes(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene
    


def Load_FPKM(Fil):
    data_type = ({'names':['gene_name' , 'FPKM'],
                   'formats':['U64' ,  np.float]})
    data = np.loadtxt(Fil , skiprows=1 , usecols=(1 , 7) , dtype = data_type)
    
    return data



F35_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F35_R1.txt')
F35_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F35_R2.txt')
F35_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F35_R3.txt')
F37_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F37_R1.txt')
F37_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F37_R2.txt')
F37_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F37_R3.txt')
F41_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F41_R1.txt')
F41_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F41_R2.txt')
F41_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\F41_R3.txt')
NT3_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT3_R1.txt')
NT3_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT3_R2.txt')
NT3_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT3_R3.txt')

NT2_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT2_R1.txt')
NT2_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT2_R2.txt')
NT2_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT2_R3.txt')
NT4_R1 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT4_R1.txt')
NT4_R2 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT4_R2.txt')
NT4_R3 = Load_FPKM('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\FPKM\\NT4_R3.txt')



union_gene = get_raw_genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
gene_name = union_gene['gene_name']

genes = []
for gene in union_gene:
    name = gene['gene_name']
    f35_1 = F35_R1[F35_R1['gene_name'] == name][0]['FPKM']
    f35_2 = F35_R2[F35_R2['gene_name'] == name][0]['FPKM']
    f35_3 = F35_R3[F35_R3['gene_name'] == name][0]['FPKM']
    f37_1 = F37_R1[F37_R1['gene_name'] == name][0]['FPKM']
    f37_2 = F37_R2[F37_R2['gene_name'] == name][0]['FPKM']
    f37_3 = F37_R3[F37_R3['gene_name'] == name][0]['FPKM']
    f41_1 = F41_R1[F41_R1['gene_name'] == name][0]['FPKM']
    f41_2 = F41_R2[F41_R2['gene_name'] == name][0]['FPKM']
    f41_3 = F41_R3[F41_R3['gene_name'] == name][0]['FPKM']
    nt2_1 = NT2_R1[NT2_R1['gene_name'] == name][0]['FPKM']
    nt2_2 = NT2_R2[NT2_R2['gene_name'] == name][0]['FPKM']
    nt2_3 = NT2_R3[NT2_R3['gene_name'] == name][0]['FPKM']
    nt3_1 = NT3_R1[NT3_R1['gene_name'] == name][0]['FPKM']
    nt3_2 = NT3_R2[NT3_R2['gene_name'] == name][0]['FPKM']
    nt3_3 = NT3_R3[NT3_R3['gene_name'] == name][0]['FPKM']
    nt4_1 = NT4_R1[NT4_R1['gene_name'] == name][0]['FPKM']
    nt4_2 = NT4_R2[NT4_R2['gene_name'] == name][0]['FPKM']
    nt4_3 = NT4_R3[NT4_R3['gene_name'] == name][0]['FPKM']
    
    genes.append((gene['gene_name'] , gene['chr'] , gene['strand'] , gene['start'] , gene['end'] , f35_1 , f35_2 , f35_3 , f37_1 , f37_2 , f37_3 , f41_1 , f41_2 , f41_3 , nt2_1 , nt2_2 , nt2_3, nt3_1 , nt3_2 , nt3_3, nt4_1 , nt4_2 , nt4_3))
    
    
genes = np.array(genes , dtype = d_type)


out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\RNA_new_gene_expression_all.txt' , 'w')
out.writelines('\t'.join(['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35_R1' , 'F35_R2' , 'F35_R3' , 'F37_R1' , 'F37_R2' , 'F37_R3' , 'F41_R1' , 'F41_R2' , 'F41_R3' , 'NT2_R1' , 'NT2_R2' , 'NT2_R3' , 'NT3_R1' , 'NT3_R2' , 'NT3_R3' , 'NT4_R1' , 'NT4_R2' , 'NT4_R3']) + '\n')

for i in gene_name:
    gene = genes[genes['gene_name'] == i][0]
    out.writelines('\t'.join([str(x) for x in gene]) + '\n')
out.close()
    


RNA_all = {'F35_R1':[] , 'F35_R2':[] , 'F35_R3':[] , 'F37_R1':[] , 'F37_R2':[] , 'F37_R3':[] ,
           'F41_R1':[] , 'F41_R2':[] , 'f41_R3':[] , 'NT3_R1':[] , 'NT3_R2':[] , 'NT3_R3':[] , 
           'NT2_R1':[] , 'NT2_R2':[] , 'NT2_R3':[] , 'NT4_R1':[] , 'NT4_R2':[] , 'NT4_R3':[]}            

RNA_all['F35_R1'] = np.log2(genes['F35_R1'] + 1)
RNA_all['F35_R2'] = np.log2(genes['F35_R2'] + 1)
RNA_all['F35_R3'] = np.log2(genes['F35_R3'] + 1)
RNA_all['F37_R1'] = np.log2(genes['F37_R1'] + 1)
RNA_all['F37_R2'] = np.log2(genes['F37_R2'] + 1)
RNA_all['F37_R3'] = np.log2(genes['F37_R3'] + 1)
RNA_all['F41_R1'] = np.log2(genes['F41_R1'] + 1)
RNA_all['F41_R2'] = np.log2(genes['F41_R2'] + 1)
RNA_all['F41_R3'] = np.log2(genes['F41_R3'] + 1)
RNA_all['NT3_R1'] = np.log2(genes['NT3_R1'] + 1)
RNA_all['NT3_R2'] = np.log2(genes['NT3_R2'] + 1)
RNA_all['NT3_R3'] = np.log2(genes['NT3_R3'] + 1)
RNA_all['NT2_R1'] = np.log2(genes['NT2_R1'] + 1)
RNA_all['NT2_R2'] = np.log2(genes['NT2_R2'] + 1)
RNA_all['NT2_R3'] = np.log2(genes['NT2_R3'] + 1)
RNA_all['NT4_R1'] = np.log2(genes['NT4_R1'] + 1)
RNA_all['NT4_R2'] = np.log2(genes['NT4_R2'] + 1)
RNA_all['NT4_R3'] = np.log2(genes['NT4_R3'] + 1)

    
    
all_cells = ['F35_R1ttF35_R2' , 'F35_R1ttF35_R3' , 'F35_R2ttF35_R3' , 'F37_R1ttF37_R2' , 'F37_R1ttF37_R3' , 'F37_R2ttF37_R3' , 
             'F41_R1ttF41_R2' , 'F41_R1ttF41_R3' , 'F41_R2ttF41_R3' , 'NT2_R1ttNT2_R2' , 'NT2_R1ttNT2_R3' , 'NT2_R2ttNT2_R3' , 
             'NT3_R1ttNT3_R2' , 'NT3_R1ttNT3_R3' , 'NT3_R2ttNT3_R3' , 'NT4_R1ttNT4_R2' , 'NT4_R1ttNT4_R3' , 'NT4_R2ttNT4_R3']
             

size = (12, 12)
Left = 0.1 ; HB = 0.1 ; width = 0.8 ; HH = 0.8

pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\Correlation_plot\\RNA_new_log2(FPKM)_all.pdf')
cor_data = {}
for c in all_cells:
    c1 = c.split("tt")[0]
    c2 = c.split("tt")[1]
    cor = round(np.corrcoef(RNA_all[c1] , RNA_all[c2])[0][1],5)
    cor_data[c] = cor
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes([Left , HB , width, HH])
    ax.scatter(RNA_all[c1] , RNA_all[c2] , alpha = 0.6 , c = 'red')
    ax.set_xlabel(c1 , size = 23)
    ax.set_ylabel(c2 , size = 23)
    ax.set_xlim(0,14)
    ax.set_ylim(0,14)
    ax.text(1 , 13 , 'R = ' + str(cor) , size = 20 )
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
    pp.savefig(fig) 
pp.close()
    
    


















