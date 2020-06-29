# -*- coding: utf-8 -*-
"""
Created on Thu Sep 05 20:11:22 2019

@author: han-luo
"""

import numpy as np
import os,re

from HiCHap.fastqPlus import Enzyme_Handle

def enzymeFind(genome, enzymeName, OutPath):
    """
        Finds restriction sites.
        Create a enzyme fragment file.Format and rule like :
        eg : 
            seq : 1(chr)  CTAGATCAATTCGATCTTAC    enzyme : GATC
            the file will be :
            1 1 4    ->   [1,4)
            1 4 13   ->   [4,13) 
            1 13 20  ->   [13,20)
        
    Parameters
    ----------
    genome : str
        Original genome file
    
    enzymeName : str
        enzyme Site
        
    """
    
    # loading Genome        
    dict_genome = {}
    file_genome = open(genome,'r')
    genome_name = genome.split('/')[-1].rstrip('.fa')
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>':
            chrs = (line.split('>')[1]).split()[0].lstrip('chr')
            dict_genome[chrs] = []
        else:
            dict_genome[chrs].extend(lists)
    # Finding enzyme site

    enzyme_site, cutsite = Enzyme_Handle(enzymeName)
    
    enzyme_file = os.path.join(OutPath,enzymeName+'_'+genome_name+'_fragments.txt')
    f = open(enzyme_file,'w')
    n = 1
    for k in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X' , 'Y']:
        seq_str = ''.join(dict_genome[k]).upper()
        pos_list = [int(m.start() + 1 + cutsite[0]) for m in re.finditer(enzyme_site, seq_str)]
        pos_list = [1] + pos_list + [len(seq_str)]
        for i in range(len(pos_list)-1):
            line = [k, str(pos_list[i]), str(pos_list[i+1])]
            f.writelines('\t'.join(line)+ '\t' + str(n) + '\n')
            n += 1
    f.close()


enzymeFind('/public/home/xxli/data/ref/mm10/mm10.fa' , 'DpnII' , './')
data_type = np.dtype({'names':['chr' , 'start' , 'end' , 'fragID'] , 
                      'formats':['S64' , np.int , np.int , np.int]})
                      
frag_data = np.loadtxt('D:\\Klf4\\cHiC\\genome\\DpnII_mm10_fragments.txt' , dtype = data_type)

frag_data_4 = frag_data[frag_data['chr'] == '4']
w = open('D:\\Klf4\\cHiC\\genome\\mm10_chr4.rmap' , 'w')
for i in frag_data_4:
    w.writelines('chr' + '\t'.join([str(x) for x in i]) + '\n')
w.close()

