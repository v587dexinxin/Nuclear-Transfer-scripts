# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 11:51:46 2017

@author: xxli
"""
import re
f1=open('E:/data/Mouse_SNP/chr1.fa','rt')
f2=open('E:/data/Mouse_SNP/DBA_2J.mgp.v5.snps.dbSNP142.vcf','rt')
f3=open('E:/data/Mouse_SNP/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf','rt')
#finalgenome=open('E:/data/Mouse_SNP/mm10.fa','wt')
ori=f1.readlines()
male_snp=f2.readlines()
fmale_snp=f3.readlines()
genome={}
male={}
fmale={}
fasta={}
#chrs=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY')
chrs=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y')
for g in chrs:
    genome=[]
    male[g]=[]
    fmale[g]=[]
for i in ori:
    if re.match('>',i):
        g=i[-2]
        fasta[g]=[]
    else:
        i=i.strip('\n')
        fasta[g].extend(list(i))
for x in male_snp:
    x=x.split()
    if male.has_key(x[0]):
        male[]
f1.close()
f2.close()
f3.close()