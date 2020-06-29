# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 19:09:22 2017

@author: DELL
"""
##   Useage: python genome_replace.py SNP.vcf GRCm38_68.fa

import numpy as np
import sys
snp=sys.argv[1]
GRCm38=sys.argv[2]
snp_type = [('chr', '<S2'), ('S', '<i4'), ('l1', '<S1'), ('l2', '<S1')]
data_snp = np.loadtxt(snp,usecols = (0,1,3,4),dtype=snp_type,skiprows=69)

def read_genome(filename):
    file_genome = open(filename)
    dict_genome = {}
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>' :
            chrs = (line.split('>')[1]).split()[0]
            dict_genome[chrs] = []
        
        else :
            dict_genome[chrs].extend(lists)
    return dict_genome

genome_m38 = read_genome(GRCm38)

new_genome = genome_m38     
for x in data_snp:
    new_genome[x[0]][int(x[1]) - 1 ] = x[3]

file_out = open('genome_C57BL.fa','w')
for k in sorted(new_genome.keys()):
    out_list = new_genome[k]
    lens = len(out_list)
    n = int(lens)/60 + 1
    file_out.write('>'+ k +' dna:chromosome chromosome:GRCm38:1:1:'+ str(lens) +':1 REF' + '\n')
    for i in range(n):
        strs = ''.join(out_list[(i*60) : ((i+1)*60)])
        file_out.write(strs+'\n')
