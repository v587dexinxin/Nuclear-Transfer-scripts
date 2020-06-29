# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:43:23 2020

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
        tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
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

chip1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw")
chip2 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349665_CC_H3K9me3_rep2.bw")
input1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw")
input2 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal_literature/GSM4349667_CC_input_rep2.bw")



Chip1 = Sig_To_50bp(chip1)
Chip2 = Sig_To_50bp(chip2)
Input1 = Sig_To_50bp(input1)
Input2 = Sig_To_50bp(input2)




pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
pc = {'c1' : [] , 'c2' : [] , 'c3' : [] , 'c4' : [] , 'c5' : [] , 'c6' : [] , 'c7' : [] , 'c8' : [] }


for c in pc:
    for g in chroms:
        tmp_pc = pc_data[c][pc_data[c]['chr'] == g]
        tmp1 = Chip1[g]
        tmp2 = Chip2[g]
        tmp3 = Input1[g]
        tmp4 = Input2[g] 
        sum1 = tmp1.sum()
        sum2 = tmp2.sum()
        sum3 = tmp3.sum()
        sum4 = tmp4.sum()
        for i in tmp_pc:
            start = i[1] // 50
            end = i[2] // 50
            data1 = tmp1[start:end].sum()  / sum1
            data2 = tmp2[start:end].sum()  / sum2
            data3 = tmp3[start:end].sum()  / sum3
            data4 = tmp4[start:end].sum()  / sum4
            if (data1 > data3) or (data2 > data4):
                pc[c].append((g , start * 50 , end * 50))
        
for c in ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']:
    print c,len(pc[c]) / len(pc_data[c])    
    out = open('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/H3K9me3_marked_Compartment_cluster'+ c[1] + '.txt' , 'w')
    for i in pc[c]:
        out.writelines('\t'.join([i[0] , str(i[1]) , str(i[2])]) + '\n')
    out.close()
    


mm_dtype = np.dtype({'names':['chr','length'] ,                                       
                     'formats':['S6',np.int]}) 
mm_10 = np.loadtxt('/public/home/xxli/data/ref/haplotype/mm10.txt' , dtype = mm_dtype)

pc = []
n = 0
for c in mm_10:
    if (c['chr'] == 'X') or (c['chr'] == 'Y'):
        continue
    g = c['chr']
    length = c['length']
    tmp1 = Chip1[g]
    tmp2 = Chip2[g]
    tmp3 = Input1[g]
    tmp4 = Input2[g] 
    sum1 = tmp1.sum()
    sum2 = tmp2.sum()
    sum3 = tmp3.sum()
    sum4 = tmp4.sum()
    for i in range(length // 200000):
        n += 1
        start = i * 200000// 50
        end = (i+1) * 200000// 50
        data1 = tmp1[start:end].sum()  / sum1
        data2 = tmp2[start:end].sum()  / sum2
        data3 = tmp3[start:end].sum()  / sum3
        data4 = tmp4[start:end].sum()  / sum4
        if (data1 > data3) or (data2 > data4):
            pc.append((g , start * 50 , end * 50))
            
            
            
            