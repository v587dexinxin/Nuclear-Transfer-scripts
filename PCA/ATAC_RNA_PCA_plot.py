# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:14:57 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import cPickle
import sys
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')



def Sig_To_100K(ATACData , mm10):
    """
    """
    
    chroms = [str(i) for i in range(1,20)] + ['X' , 'Y']
    New_Data = np.array([])
    for c in chroms:
        tmp_data = ATACData[ATACData['chr'] == c]
        max_ = mm10[mm10['chr'] == c][0]['length']
        bin_size = max_ // 100000 + 1
        tmp = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100000
            end = line['end'] // 100000
            tmp[start] += line['score']
            for i in range(start + 1 , end):
                tmp[i] += line['score']
        New_Data = np.hstack((New_Data , tmp))
    New_Data = np.log(New_Data / 100 + 1)
    return New_Data


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
#ATAC
lenth_type = np.dtype({'names':['chr' , 'length'],
                      'formats':['S64' , np.int ]})                      
mm10 = np.loadtxt('/public/home/xxli/data/ref/haplotype/mm10.txt' , dtype = lenth_type)

sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                  'formats':['S4',np.int , np.int , np.float]})

CCSData = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/CCS_ATAC_1K.bedgraph' , dtype = sig_type)
NT5Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT5_ATAC_1K.bedgraph' , dtype = sig_type)
NT6Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT6_ATAC_1K.bedgraph' , dtype = sig_type)
fESCData = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/fESC_ATAC_1K.bedgraph' , dtype = sig_type)

CCS_R1Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/CCS_R1_ATAC_1K.bedgraph', dtype = sig_type)
CCS_R2Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/CCS_R2_ATAC_1K.bedgraph', dtype = sig_type)
NT5_R1Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT5_R1_ATAC_1K.bedgraph', dtype = sig_type)
NT5_R2Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT5_R2_ATAC_1K.bedgraph', dtype = sig_type)
NT6_R1Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT6_R1_ATAC_1K.bedgraph', dtype = sig_type)
NT6_R2Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/NT6_R2_ATAC_1K.bedgraph', dtype = sig_type)
fESC_R1Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/fESC_R1_ATAC_1K.bedgraph', dtype = sig_type)
fESC_R2Data = np.loadtxt('/public/home/xxli/data/BDF1_New/ATAC/workspace/signal/normalization/bedgraph_1K/fESC_R2_ATAC_1K.bedgraph', dtype = sig_type)


CCS_R1_100K = Sig_To_100K(CCS_R1Data , mm10)
CCS_R2_100K = Sig_To_100K(CCS_R2Data , mm10)
NT5_R1_100K = Sig_To_100K(NT5_R1Data , mm10)
NT5_R2_100K = Sig_To_100K(NT5_R2Data , mm10)
NT6_R1_100K = Sig_To_100K(NT6_R1Data , mm10)
NT6_R2_100K = Sig_To_100K(NT6_R2Data , mm10)
fESC_R1_100K = Sig_To_100K(fESC_R1Data , mm10)
fESC_R2_100K = Sig_To_100K(fESC_R2Data , mm10)

CCS_100K = Sig_To_100K(CCSData , mm10)
NT5_100K = Sig_To_100K(NT5Data , mm10)
NT6_100K = Sig_To_100K(NT6Data , mm10)
fESC_100K = Sig_To_100K(fESCData , mm10)




matrix = np.vstack((CCS_R1_100K , CCS_R2_100K , NT5_R1_100K , NT5_R2_100K , NT6_R1_100K , NT6_R2_100K , fESC_R1_100K , fESC_R2_100K))
matrix = np.vstack((CCS_100K , NT5_100K , NT6_100K , fESC_100K))
pca = PCA(n_components=2)
pca.fit(matrix)
reco = np.matrix(matrix)*np.matrix(pca.components_.T)
print(pca.explained_variance_ratio_)

#Results
'''
Replicates
print(pca.explained_variance_ratio_) : [0.69793182 0.16611366]
matrix([[ 11.64806898, -14.82449573],
        [  9.49057371, -11.66542195],
        [ -4.77163177, -15.75784752],
        [ -4.36699655,  -6.11966147],
        [ -3.28221468, -13.55243737],
        [ -5.07352968, -11.10967604],
        [ -5.70941931, -15.13150442],
        [ -5.38651679, -17.13828416]])
        
Merge
print(pca.explained_variance_ratio_) : [0.87261022 0.10243619]
matrix([[11.03726062,  9.51686558],
        [-4.22907534,  6.63076656],
        [-3.943351  ,  8.11338475],
        [-5.15572736, 12.8516227 ]])
'''
pc1 = 0.69793182 ; pc2 = 0.16611366
reco = np.matrix([[ 11.64806898, -14.82449573],
        [  9.49057371, -11.66542195],
        [ -4.77163177, -15.75784752],
        [ -4.36699655,  -6.11966147],
        [ -3.28221468, -13.55243737],
        [ -5.07352968, -11.10967604],
        [ -5.70941931, -15.13150442],
        [ -5.38651679, -17.13828416]])
x = np.array(reco.T[0])[0]
y = np.array(reco.T[1])[0]
left, bottom, width, height = 0.2, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.scatter(x[:2] , y[:2] , color = 'darkorange' , s = 300)
ax.text(x[:2].min() - 0.1 , y[:2].min() ,'CCS')

ax.scatter(x[2:4] , y[2:4] , color = 'red' , s =300)
ax.text(x[2:4].min() + 0.1 , y[2:4].min() ,'NT5')

ax.scatter(x[4:6] , y[4:6] , color = 'green' , s = 300)
ax.text(x[4:6].min() + 0.1 , y[4:6].min() ,'NT6')

ax.scatter(x[6:] , y[6:] , color = 'blue' , s = 300)
ax.text(x[6:].min() + 0.1 , y[6:].min() ,'fESC')

ax.set_xlim((-7 , 12.5))
ax.set_ylim((-18 , -5))
ax.set_xlabel('PC1:' + str(pc1))
ax.set_ylabel('PC2:' + str(pc2))
ax.set_title('ATAC')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\ATAC_PCA_replicates.pdf')


pc1 = 0.87261022 ; pc2 = 0.10243619
reco = np.matrix([[11.03726062,  9.51686558],
        [-4.22907534,  6.63076656],
        [-3.943351  ,  8.11338475],
        [-5.15572736, 12.8516227 ]])
x = np.array(reco.T[0])[0]
y = np.array(reco.T[1])[0]
left, bottom, width, height = 0.2, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.scatter(x[0] , y[0] , color = 'darkorange' , s = 300)
ax.text(x[0].min() - 0.1 , y[0].min() ,'CCS')

ax.scatter(x[1] , y[1] , color = 'red' , s = 300)
ax.text(x[1].min() + 0.1 , y[1].min() ,'NT5')

ax.scatter(x[2] , y[2] , color = 'green' , s = 300)
ax.text(x[2].min() + 0.1 , y[2].min() ,'NT6')

ax.scatter(x[3] , y[3] , color = 'blue' , s = 300)
ax.text(x[3].min() + 0.1 , y[3].min() ,'fESC')
ax.set_xlabel('PC1:' + str(pc1))
ax.set_ylabel('PC2:' + str(pc2))
ax.set_title('ATAC')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\ATAC_PCA_Merged.pdf')

#RNA
union_type = np.dtype({'names':['gene_name' , 'chr' , 'starnd' , 'start' , 'end' , 'CCS_R1' , 'CCS_R2' , 'CCS_R3' , 'NT5_R1' , 'NT5_R2' , 'NT5_R3' , 'NT5_R4' , 'NT6_R1' , 'NT6_R2' , 'NT6_R3' , 'fESC_R1' , 'fESC_R2' , 'fESC_R3'],
                       'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_gene = np.loadtxt('/public/home/xxli/data/BDF1_New/RNA/workspace/gene_expression/all_genes.txt' , skiprows = 1 , dtype = union_type , usecols = (1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 9 , 10 , 11 ))

#All gene
matrix = []
for i in union_gene:
    if (i['CCS_R1'] > 0) or (i['CCS_R2'] > 0) or (i['CCS_R3'] > 0) or (i['NT5_R1'] > 0) or (i['NT5_R2'] > 0) or (i['NT5_R3'] > 0) or (i['NT5_R4'] > 0) or (i['NT6_R1'] > 0) or (i['NT6_R2'] > 0) or (i['NT6_R3'] > 0) or (i['fESC_R1'] > 0) or (i['fESC_R2'] > 0) or (i['fESC_R3'] > 0):
        matrix.append((i['CCS_R1'] , i['CCS_R2'] , i['CCS_R3'] , i['NT5_R1'] , i['NT5_R2'] , i['NT5_R3'] , i['NT5_R4'] , i['NT6_R1'] , i['NT6_R2'] , i['NT6_R3'] , i['fESC_R1'] , i['fESC_R2'] , i['fESC_R3']))

'''
All genes
print(pca.explained_variance_ratio_) : [0.74762914 0.10625292]
matrix([[ 56.36893769,  19.74780811],
        [ 52.78296881,  16.67754091],
        [ 56.02230224,  17.57937572],
        [-78.50364173,  31.90355237],
        [-68.50719233,  53.23124153],
        [-60.99500619,  60.61930666],
        [-78.84859007,  24.21018893],
        [-78.56593597,  29.62743462],
        [-81.0896022 ,  24.93942114],
        [-80.0979529 ,  -9.44137101],
        [-80.72317242, -10.83165628],
        [-80.24860035,  -6.96357292]])
'''
pc1 = 0.74762914 ; pc2 = 0.10625292
reco =  np.matrix([[ 56.36893769,  19.74780811],
        [ 52.78296881,  16.67754091],
        [ 56.02230224,  17.57937572],
        [-78.50364173,  31.90355237],
        [-68.50719233,  53.23124153],
        [-60.99500619,  60.61930666],
        [-78.84859007,  24.21018893],
        [-78.56593597,  29.62743462],
        [-81.0896022 ,  24.93942114],
        [-80.0979529 ,  -9.44137101],
        [-80.72317242, -10.83165628],
        [-80.24860035,  -6.96357292]])
        
x = np.array(reco.T[0])[0]
y = np.array(reco.T[1])[0]
left, bottom, width, height = 0.2, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.scatter(x[:3] , y[:3] , color = 'darkorange' , s = 300)
ax.text(x[:3].min() - 10 , y[:3].min() ,'CCS')

ax.scatter(x[3:6] , y[3:6] , color = 'red' , s = 300)
ax.text(x[3:6].min() + 10 , y[3:6].min() ,'NT5')

ax.scatter(x[6:9] , y[6:9] , color = 'green' , s = 300)
ax.text(x[6:9].min() + 10 , y[6:9].min() ,'NT6')

ax.scatter(x[9:] , y[9:] , color = 'blue' , s = 300)
ax.text(x[9:].min() + 10 , y[9:].min() ,'fESC')

ax.set_xlim((-90 , 70))
ax.set_ylim((-15 , 65))
ax.set_xlabel('PC1:' + str(pc1))
ax.set_ylabel('PC2:' + str(pc2))
ax.set_title('Gene expression')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\PCA_Gene_expression_Allgenes_replicates.pdf')



        
#diff expression        
matrix = []
for i in union_gene:
    if (i['CCS_R1'] > 0) or (i['CCS_R2'] > 0) or (i['CCS_R3'] > 0) or (i['NT5_R1'] > 0) or (i['NT5_R3'] > 0) or (i['NT5_R4'] > 0) or (i['NT6_R1'] > 0) or (i['NT6_R2'] > 0) or (i['NT6_R3'] > 0) or (i['fESC_R1'] > 0) or (i['fESC_R2'] > 0) or (i['fESC_R3'] > 0):
        CCS_F = (i['CCS_R1'] + i['CCS_R2'] + i['CCS_R3']) / 3
        NT5_F = (i['NT5_R1'] + i['NT5_R3'] + i['NT5_R4']) / 3
        NT6_F = (i['NT6_R1'] + i['NT6_R2'] + i['NT6_R3']) / 3
        fESC_F = (i['fESC_R1'] + i['fESC_R2'] + i['fESC_R3']) / 3
        if ((CCS_F > 1.5 * NT5_F) and (CCS_F > 1.5 * NT6_F) and (CCS_F > 1.5 * fESC_F)) or ((1.5 * CCS_F < NT5_F) and (1.5 * CCS_F < NT6_F) and (1.5 * CCS_F < fESC_F)):
            matrix.append([i['CCS_R1'] , i['CCS_R2'] , i['CCS_R3'] , i['NT5_R1'] , i['NT5_R3'] , i['NT5_R4'] , i['NT6_R1'] , i['NT6_R2'] , i['NT6_R3'] , i['fESC_R1'] , i['fESC_R2'] , i['fESC_R3']])
matrix = np.array(matrix).T
matrix = np.log2(matrix + 1)
        
pca = PCA(n_components=2)
pca.fit(matrix)
reco = np.matrix(matrix)*np.matrix(pca.components_.T)
print(pca.explained_variance_ratio_)
        
#Results
'''       
Diff expression
print(pca.explained_variance_ratio_) : [0.8442915 0.0671485]
matrix([[  78.24578218,   27.85186704],
        [  73.34981765,   22.8816493 ],
        [  78.39688432,   22.77347094],
        [-111.92636949,   36.88399958],
        [ -98.70407019,   62.99043649],
        [ -88.62916614,   70.6227921 ],
        [-112.19729489,   27.64016256],
        [-112.30895161,   36.98769168],
        [-115.7431809 ,   30.55719811],
        [-113.92428326,   -2.44506425],
        [-113.89989952,   -5.89506471],
        [-113.7841269 ,   -1.85234979]])

''' 
pc1 =    0.8442915 ; pc2 = 0.0671485
reco =  np.matrix([[  78.24578218,   27.85186704],
        [  73.34981765,   22.8816493 ],
        [  78.39688432,   22.77347094],
        [-111.92636949,   36.88399958],
        [ -98.70407019,   62.99043649],
        [ -88.62916614,   70.6227921 ],
        [-112.19729489,   27.64016256],
        [-112.30895161,   36.98769168],
        [-115.7431809 ,   30.55719811],
        [-113.92428326,   -2.44506425],
        [-113.89989952,   -5.89506471],
        [-113.7841269 ,   -1.85234979]])
        
x = np.array(reco.T[0])[0]
y = np.array(reco.T[1])[0]
left, bottom, width, height = 0.2, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.scatter(x[:3] , y[:3] , color = 'darkorange' , s = 70)
ax.text(x[:3].min() - 10 , y[:3].min() ,'CCS')

ax.scatter(x[3:6] , y[3:6] , color = 'red' , s = 70)
ax.text(x[3:6].min() + 10 , y[3:6].min() ,'NT5')

ax.scatter(x[6:9] , y[6:9] , color = 'green' , s = 70)
ax.text(x[6:9].min() + 10 , y[6:9].min() ,'NT6')

ax.scatter(x[9:] , y[9:] , color = 'blue' , s = 70)
ax.text(x[9:].min() + 10 , y[9:].min() ,'fESC')

#ax.set_xlim((-90 , 70))
#ax.set_ylim((-15 , 65))
ax.set_xlabel('PC1:' + str(pc1))
ax.set_ylabel('PC2:' + str(pc2))
ax.set_title('Gene expression')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\PCA_Gene_expression_Diffgenes_replicates.pdf')




#Merged
matrix = []
for i in union_gene:
    if (i['CCS_R1'] > 0) or (i['CCS_R2'] > 0) or (i['CCS_R3'] > 0) or (i['NT5_R1'] > 0) or (i['NT5_R2'] > 0) or (i['NT5_R3'] > 0) or (i['NT5_R4'] > 0) or (i['NT6_R1'] > 0) or (i['NT6_R2'] > 0) or (i['NT6_R3'] > 0) or (i['fESC_R1'] > 0) or (i['fESC_R2'] > 0) or (i['fESC_R3'] > 0):
        CCS_F = (i['CCS_R1'] + i['CCS_R2'] + i['CCS_R3']) / 3
        NT5_F = (i['NT5_R1'] + i['NT5_R3'] + i['NT5_R4']) / 3
        NT6_F = (i['NT6_R1'] + i['NT6_R2'] + i['NT6_R3']) / 3
        fESC_F = (i['fESC_R1'] + i['fESC_R2'] + i['fESC_R3']) / 3
        matrix.append((CCS_F , NT5_F , NT6_F , fESC_F))
        
matrix = np.array(matrix).T
matrix = np.log2(matrix + 1)
        
pca = PCA(n_components=2)
pca.fit(matrix)
reco = np.matrix(matrix)*np.matrix(pca.components_.T)
print(pca.explained_variance_ratio_)

'''
print(pca.explained_variance_ratio_) : [0.82680764 0.12205374]
matrix([[  80.90650944,  -32.17231559],
        [ -99.59509573,  -77.52108686],
        [-114.02306967,  -43.90756572],
        [-115.2230636 ,   11.33545231]])
        
'''
pc1 = 0.82680764 ; pc2 = 0.12205374
reco = np.matrix([[  80.90650944,  -32.17231559],
        [ -99.59509573,  -77.52108686],
        [-114.02306967,  -43.90756572],
        [-115.2230636 ,   11.33545231]])
        
x = np.array(reco.T[0])[0]
y = np.array(reco.T[1])[0]
left, bottom, width, height = 0.2, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.scatter(x[0] , y[0] , color = 'darkorange' , s = 300)
ax.text(x[0].min() - 0.1 , y[0].min() ,'CCS')

ax.scatter(x[1] , y[1] , color = 'red' , s = 300)
ax.text(x[1].min() + 0.1 , y[1].min() ,'NT5')

ax.scatter(x[2] , y[2] , color = 'green' , s = 300)
ax.text(x[2].min() + 0.1 , y[2].min() ,'NT6')

ax.scatter(x[3] , y[3] , color = 'blue' , s = 300)
ax.text(x[3].min() + 0.1 , y[3].min() ,'fESC')
ax.set_xlabel('PC1:' + str(pc1))
ax.set_ylabel('PC2:' + str(pc2))
ax.set_title('Gene expression')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\PCA_Gene_expression_Allgenes_Merged.pdf')




