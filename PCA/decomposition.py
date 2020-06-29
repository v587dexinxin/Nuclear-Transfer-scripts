# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 20:43:54 2017

@author: DELL
"""
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
def read_file(files):
    filein = open(files)
    l = []
    next(filein)
    for line in filein:
        line = line.strip('\n')
        lists = line.split('\t')
        l.append(np.log(float(lists[9])/1000+1))
    return np.array(l)

CCS = read_file('CCS_ATAC_summary.bed')    
ESC = read_file('ESC_ATAC_summary.bed')
NT5 = read_file('NT5_ATAC_summary.bed')
NT6 = read_file('NT6_ATAC_summary.bed')
CCS_r = read_file('CCS_RNA_summary.bed')    
ESC_r = read_file('ESC_RNA_summary.bed')
NT5_r = read_file('NT5_RNA_summary.bed')
NT6_r = read_file('NT6_RNA_summary.bed')

X = np.array([CCS,ESC,NT5,NT6])
X_RNA = np.array([CCS_r,ESC_r,NT5_r,NT6_r])
pca_atac = PCA(n_components=4)
pca_rna = PCA(n_components=4)
pca_atac.fit(X)
pca_rna.fit(X_RNA)
reco_rna = np.matrix(X_RNA)*np.matrix(pca_rna.components_.T)
reco_atac = np.matrix(X)*np.matrix(pca_atac.components_.T)

plt.plot(reco_rna[:,:1],reco_rna[:,1:2],'o')
plt.text(reco_rna[0,0]-25,reco_rna[0,1],'CCS')
plt.text(reco_rna[1,0]+10,reco_rna[1,1],'ESC')
plt.text(reco_rna[2,0]+10,reco_rna[2,1],'NT5')
plt.text(reco_rna[3,0]+10,reco_rna[3,1],'NT6')
plt.title('RNA_seq')
plt.xlabel('PC1')
plt.ylabel('PC2')

plt.plot(reco_atac[:,:1],reco_atac[:,1:2],'o')
plt.text(reco_atac[0,0]-15,reco_atac[0,1],'CCS')
plt.text(reco_atac[1,0]+10,reco_atac[1,1],'ESC')
plt.text(reco_atac[2,0]+10,reco_atac[2,1],'NT5')
plt.text(reco_atac[3,0]+10,reco_atac[3,1]-5,'NT6')
plt.title('ATAC_seq')
plt.xlabel('PC1')
plt.ylabel('PC2')