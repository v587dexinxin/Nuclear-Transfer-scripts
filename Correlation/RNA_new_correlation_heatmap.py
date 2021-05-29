# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 14:45:47 2020

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
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Our Own Color Map
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')

d_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35_R1' , 'F35_R2' , 'F35_R3' , 'F37_R1' , 'F37_R2' , 'F37_R3' , 'F41_R1' , 'F41_R2' , 'F41_R3' , 'NT2_R1' , 'NT2_R2' , 'NT2_R3' , 'NT3_R1' , 'NT3_R2' , 'NT3_R3' , 'NT4_R1' , 'NT4_R2' , 'NT4_R3'],
                   'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})



                

def get_raw_genes(Fil):
    '''
    '''
    
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene
    

def get_RNA_new_genes(Fil): 
    '''
    '''
    d_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35_R1' , 'F35_R2' , 'F35_R3' , 'F37_R1' , 'F37_R2' , 'F37_R3' , 'F41_R1' , 'F41_R2' , 'F41_R3' , 'NT2_R1' , 'NT2_R2' , 'NT2_R3' , 'NT3_R1' , 'NT3_R2' , 'NT3_R3' , 'NT4_R1' , 'NT4_R2' , 'NT4_R3'],
                   'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


    union_gene = np.loadtxt(Fil , dtype = d_type , skiprows = 1)

    return union_gene


def get_ave_FPKM(union_gene):
    '''
    '''
    d_type_1 = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'F35' , 'F37' , 'F41' , 'NT2' , 'NT3' , 'NT4'],
                         'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float]})

    genes_1 = []
    for i in range(len(union_gene)):
        f35 = (union_gene[i]['F35_R1'] + union_gene[i]['F35_R2'] + union_gene[i]['F35_R3']) / 3
        f37 = (union_gene[i]['F37_R1'] + union_gene[i]['F37_R2'] + union_gene[i]['F37_R3']) / 3
        f41 = (union_gene[i]['F41_R1'] + union_gene[i]['F41_R2'] + union_gene[i]['F41_R3']) / 3
        nt2 = (union_gene[i]['NT2_R1'] + union_gene[i]['NT2_R2'] + union_gene[i]['NT2_R3']) / 3
        nt3 = (union_gene[i]['NT3_R1'] + union_gene[i]['NT3_R2'] + union_gene[i]['NT3_R3']) / 3
        nt4 = (union_gene[i]['NT4_R1'] + union_gene[i]['NT4_R2'] + union_gene[i]['NT4_R3']) / 3 
        
        genes_1.append((union_gene[i]['gene_name'] , union_gene[i]['chr'] , union_gene[i]['strand'] , 
                        union_gene[i]['start'] , union_gene[i]['end'] , f35 , f37 , f41 , nt2 , nt3 , nt4))
    genes_1 = np.array(genes_1 , dtype = d_type_1)
    
    return genes_1
    


union_gene = get_raw_genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union_gene_new = get_RNA_new_genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_new\\FPKM\\RNA_new_gene_expression_all.txt')
union_gene_new_ave = get_ave_FPKM(union_gene_new)
 




cells = {'CCS':0 , 'NT3':1 , 'NT2':2 , 'NT4':3 , 'NT5':4 , 'NT6':5 , 'F40':6 , 'F35':7 , 'F37':8 , 'F41':9}
cell = ['CCS' , 'NT3' , 'NT2' , 'NT4' , 'NT5' , 'NT6' , 'F40' , 'F35' , 'F37' , 'F41']

gene_name = union_gene['gene_name']

RNA_all = {cell[0]:[] , cell[1]:[] , cell[2]:[] , cell[3]:[] , cell[4]:[] , cell[5]:[] , cell[6]:[] , cell[7]:[] , cell[8]:[] , cell[9]:[]}

for i in range(len(gene_name)):
    if (union_gene[i]['CCS'] != 0) or (union_gene[i]['NT5'] != 0) or (union_gene[i]['NT6'] != 0) or (union_gene[i]['fESC'] != 0) or (union_gene_new_ave[i]['F35'] != 0) or (union_gene_new_ave[i]['F37'] != 0) or (union_gene_new_ave[i]['F41'] != 0) or (union_gene_new_ave[i]['NT3'] != 0) or (union_gene_new_ave[i]['NT2'] != 0) or (union_gene_new_ave[i]['NT4'] != 0):
        if (union_gene[i]['CCS'] <= 3000) and (union_gene[i]['NT5'] <= 3000) and (union_gene[i]['NT6'] <= 3000) and (union_gene[i]['fESC'] <= 3000) and (union_gene_new_ave[i]['F35'] <= 3000) and (union_gene_new_ave[i]['F37'] <= 3000) and (union_gene_new_ave[i]['F41'] <= 3000) and (union_gene_new_ave[i]['NT3'] <= 3000) and (union_gene_new_ave[i]['NT2'] <= 3000) and (union_gene_new_ave[i]['NT4'] <= 3000):
            RNA_all['CCS'].append(union_gene[i]['CCS'])
            RNA_all['NT5'].append(union_gene[i]['NT5'])
            RNA_all['NT6'].append(union_gene[i]['NT6'])
            RNA_all['F40'].append(union_gene[i]['fESC'])
            RNA_all['F35'].append(union_gene_new_ave[i]['F35'])
            RNA_all['F37'].append(union_gene_new_ave[i]['F37'])
            RNA_all['F41'].append(union_gene_new_ave[i]['F41'])
            RNA_all['NT3'].append(union_gene_new_ave[i]['NT3'])
            RNA_all['NT2'].append(union_gene_new_ave[i]['NT2'])
            RNA_all['NT4'].append(union_gene_new_ave[i]['NT4'])
            
            
for k , v in RNA_all.items():
    v = np.array(v)
    RNA_all[k] = v            


#---------------------------------Get_correlation_matrix-------------------------------------------
    
cor_matrix = np.zeros((len(cell),len(cell)))
for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(RNA_all[i] , RNA_all[j])[0][1] , 3)


pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\Correlation_plot\\RNA_new_diff_cellline_all.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=1,vmin = 0.65,cmap=my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)
for i in range(len(x)):
    ax.plot([-0.5 , len(x) - 0.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , len(x) - 0.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 14)


ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()




##hierarchy_cluster
data = np.vstack((RNA_all['CCS'] , RNA_all['NT3'] , RNA_all['NT2'] , RNA_all['NT4'] , RNA_all['NT5'] , RNA_all['NT6'] , RNA_all['F40'] , RNA_all['F35'] , RNA_all['F37'] , RNA_all['F41']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)








##--------------------------------------------------------------------------------
cells = {'CCS':0 , 'NT3':4 , 'NT2':6 , 'NT4':5 , 'NT5':2 , 'NT6':3 , 'F40':1 , 'F35':8 , 'F37':7 , 'F41':9}
cell = ['CCS' , 'NT3' , 'NT2' , 'NT4' , 'NT5' , 'NT6' , 'F40' , 'F35' , 'F37' , 'F41']
cell = ['CCS' , 'F40' , 'NT5' , 'NT6' , 'NT3' , 'NT4' , 'NT2' , 'F37' , 'F35' , 'F41']

gene_name = union_gene['gene_name']

RNA_all = {cell[0]:[] , cell[1]:[] , cell[2]:[] , cell[3]:[] , cell[4]:[] , cell[5]:[] , cell[6]:[] , cell[7]:[] , cell[8]:[] , cell[9]:[]}

for i in range(len(gene_name)):
    if (union_gene[i]['CCS'] != 0) or (union_gene[i]['NT5'] != 0) or (union_gene[i]['NT6'] != 0) or (union_gene[i]['fESC'] != 0) or (union_gene_new_ave[i]['F35'] != 0) or (union_gene_new_ave[i]['F37'] != 0) or (union_gene_new_ave[i]['F41'] != 0) or (union_gene_new_ave[i]['NT3'] != 0) or (union_gene_new_ave[i]['NT2'] != 0) or (union_gene_new_ave[i]['NT4'] != 0):
        if (union_gene[i]['CCS'] <= 3000) and (union_gene[i]['NT5'] <= 3000) and (union_gene[i]['NT6'] <= 3000) and (union_gene[i]['fESC'] <= 3000) and (union_gene_new_ave[i]['F35'] <= 3000) and (union_gene_new_ave[i]['F37'] <= 3000) and (union_gene_new_ave[i]['F41'] <= 3000) and (union_gene_new_ave[i]['NT3'] <= 3000) and (union_gene_new_ave[i]['NT2'] <= 3000) and (union_gene_new_ave[i]['NT4'] <= 3000):
            RNA_all['CCS'].append(union_gene[i]['CCS'])
            RNA_all['NT5'].append(union_gene[i]['NT5'])
            RNA_all['NT6'].append(union_gene[i]['NT6'])
            RNA_all['F40'].append(union_gene[i]['fESC'])
            RNA_all['F35'].append(union_gene_new_ave[i]['F35'])
            RNA_all['F37'].append(union_gene_new_ave[i]['F37'])
            RNA_all['F41'].append(union_gene_new_ave[i]['F41'])
            RNA_all['NT3'].append(union_gene_new_ave[i]['NT3'])
            RNA_all['NT2'].append(union_gene_new_ave[i]['NT2'])
            RNA_all['NT4'].append(union_gene_new_ave[i]['NT4'])
            
            
for k , v in RNA_all.items():
    v = np.array(v)
    RNA_all[k] = v            


#---------------------------------Get_correlation_matrix-------------------------------------------
    
cor_matrix = np.zeros((len(cell),len(cell)))
for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(RNA_all[i] , RNA_all[j])[0][1] , 3)


pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA_new\\Correlation_plot\\RNA_new_diff_cellline_all_1.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=1,vmin = 0.65,cmap=my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)
for i in range(len(x)):
    ax.plot([-0.5 , len(x) - 0.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , len(x) - 0.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 14)


ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()


['CCS' , 'NT3' , 'NT2' , 'NT4' , 'NT5' , 'NT6' , 'F40' , 'F35' , 'F37' , 'F41']

##hierarchy_cluster
data = np.vstack((RNA_all['CCS'] , RNA_all['NT3'] , RNA_all['NT2'] , RNA_all['NT4'] , RNA_all['NT5'] , RNA_all['NT6'] , RNA_all['F40'] , RNA_all['F35'] , RNA_all['F37'] , RNA_all['F41']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)







##---------------------------------Get_remove_batch_effect_matrix-----------------------

##使用limma的removeBatchEffect去除批次效应
cell = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
cor_matrix = np.array([[1.0000000 , 0.8692198 , 0.8585822 , 0.9225729 , 0.8552458] , 
[0.8692198 , 1.0000000 , 0.9939547 , 0.9899698 , 0.9836088] , 
[0.8585822 , 0.9939547 , 1.0000000 , 0.9881710 , 0.9873628] ,
[0.9225729 , 0.9899698 , 0.9881710 , 1.0000000 , 0.9846348] ,
[0.8552458 , 0.9836088 , 0.9873628 , 0.9846348 , 1.0000000]])

cor_matrix = np.round(cor_matrix , 3)

pp = PdfPages('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\Correlation_plot\\RNA_correlation_heatmap_remove_batch_1.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=0.98,vmin = 0.91,cmap=my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)
for i in range(len(x)):
    ax.plot([-0.5 , len(x) - 0.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , len(x) - 0.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 14)


ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()















