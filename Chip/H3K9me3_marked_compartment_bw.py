# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:59:07 2020

@author: han-luo
"""

from __future__ import division
import numpy as np 
import pyBigWig



def Sig_To_50bp(signal):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals(g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 50 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 50
            end = line['end'] // 50
            for i in range(start,end):
                New_Data[g][i] += line['value']
    
    return New_Data
    
    
    
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


                    
pc1 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc2 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc4 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc5 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc6 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc7 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc8 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )

chip1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_chip_50bp.bw")
input1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_input_50bp.bw")


Chip1 = Sig_To_50bp(chip1)
Input1 = Sig_To_50bp(input1)




pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
pc = {'c1' : [] , 'c2' : [] , 'c3' : [] , 'c4' : [] , 'c5' : [] , 'c6' : [] , 'c7' : [] , 'c8' : [] }


for c in pc:
    for g in chroms:
        tmp_pc = pc_data[c][pc_data[c]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_pc:
            start = i[1] // 50
            end = i[2] // 50
            data1 = tmp1[start:end].sum() / sum1
            data3 = tmp3[start:end].sum()  / sum3
            if data1 > data3:
                pc[c].append((g , start * 50 , end * 50))
        
for c in ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']:
    print c,len(pc[c]) / len(pc_data[c])   
