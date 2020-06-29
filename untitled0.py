# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 17:23:02 2020

@author: han-luo
"""

import numpy as np 
import pyBigWig


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

sum1 = chip1.header()['sumData'] 
sum2 = chip2.header()['sumData']
sum3 = input1.header()['sumData'] 
sum4 = input2.header()['sumData']

pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
pc = {'c1' : [] , 'c2' : [] , 'c3' : [] , 'c4' : [] , 'c5' : [] , 'c6' : [] , 'c7' : [] , 'c8' : [] }


for c in pc:
    for g in chroms:
        tmp_pc = pc_data[c][pc_data[c]['chr'] == g]
        tmp1 = np.array(list(chip1.intervals('chr' + g)) , dtype = signal_type)
        tmp2 = np.array(list(chip2.intervals('chr' + g)) , dtype = signal_type)
        tmp3 = np.array(list(input1.intervals('chr' + g)) , dtype = signal_type)
        tmp4 = np.array(list(input2.intervals('chr' + g)) , dtype = signal_type)    
        for i in tmp_pc:
            start = i[1]
            end = i[2]
            mask1 = (tmp1['start'] < end) & (tmp1['end'] > start)
            mask2 = (tmp2['start'] < end) & (tmp2['end'] > start)
            mask3 = (tmp3['start'] < end) & (tmp3['end'] > start)
            mask4 = (tmp4['start'] < end) & (tmp4['end'] > start)
            overlap1 = tmp1[mask1]
            overlap2 = tmp2[mask2]
            overlap3 = tmp3[mask3]
            overlap4 = tmp4[mask4]
            data1 = overlap1['value'].sum() / sum1
            data2 = overlap2['value'].sum() / sum2
            data3 = overlap3['value'].sum() / sum3
            data4 = overlap4['value'].sum() / sum4
            if (data1 > data3) or (data2 > data4):
                pc[c].append((g , start , end))
        
        
        
for c in pc:
    print c,len(pc[c]) / len(pc_data[c])
    
    
    
    
    
    
    
    