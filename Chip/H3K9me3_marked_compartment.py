# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 15:07:33 2020

@author: xxli
"""


from __future__ import division
import numpy as np 
import pyBigWig
import pandas as pd
import matplotlib 
# Use a non-interactive backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import islice  
    
    
def Sig_To_100bp(signal , control):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        if control == 'input':
            tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        else:
            tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            # end = line['end'] // 100
            New_Data[g][start] += line['value']
    
    return New_Data
    


def get_union_gene_sites(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS_FPKM'] , i['NT2_FPKM'] , i['NT3_FPKM'] , i['NT4_FPKM'] , i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F37_FPKM'] , i['F40_FPKM'] , i['F41_FPKM']))
    union = np.array(union , dtype = union_type_1)
    return union


def get_raw_genes_new(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT2' , 'NT3' , 'NT4' , 'NT5' , 'NT6' , 'F35' , 'F37' , 'F40' , 'F41'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene
    


def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf


def Get_compartment_genes(pc , union_gene):
    '''
    '''
    genes = []
    for i in pc:
        g = i[0]
        start = i[1] 
        end = i[2]
        union_tmp = union_gene[union_gene['chr'] == 'chr' + g]
        mask = (union_tmp['gene_site'] >= start) & (union_tmp['gene_site'] <= end)
        overlap = union_tmp[mask]
        if overlap.size != 0:
            for j in overlap:
                genes.append(j)
            
            
    genes = np.array(genes , dtype = union_gene.dtype)
    return genes


def OutputTXT(Info , outfil , head):
    """
    """
    print ("Output ...")
    
    # head = ['chr','position','PC-M','PC-P','diff','P_Value']
    with open(outfil, 'w') as out:
        out.writelines('\t'.join(head)+'\n')
        
        for line in Info:
            line = map(str, line)
            out.writelines('\t'.join(line)+'\n')
    out.close()
    
                
def Write2fils_nochr(filname , peaks):    
    '''
    '''
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
    
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})
pc_type_1 = np.dtype({'names':['chr' , 'pc'],
                      'formats':['S8' , np.float]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000





                    
pc1 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc2 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc4 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc5 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc6 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc7 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc8 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
ccs = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/CCS_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type_1)
ccs = ccs[ccs['chr'] != 'X']

# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")

# chip1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_chip_50bp.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal/uniq_pairs_CCS_H3K9me3_input_50bp.bw")


chip1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw")
input1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw")


Chip1 = Sig_To_100bp(chip1 , 'chip')
Input1 = Sig_To_100bp(input1 , 'input')




pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
pc = {'c1' : [] , 'c2' : [] , 'c3' : [] , 'c4' : [] , 'c5' : [] , 'c6' : [] , 'c7' : [] , 'c8' : [] }


n = 0
for g in chroms:
    tmp = ccs[ccs['chr'] == g]
    tmp1 = Chip1[g]
    tmp3 = Input1[g]
    sum1 = tmp1.sum()
    sum3 = tmp3.sum()
    for i in range(len(tmp)):
        start = i * 200000 // 100
        end = (i + 1) * 200000 // 100
        data1 = tmp1[start:end].sum() 
        data3 = tmp3[start:end].sum()
        if data1 > 1.4 * data3:
            n += 1

for c in pc:
    for g in chroms:
        tmp_pc = pc_data[c][pc_data[c]['chr'] == g]
        tmp1 = Chip1[g]
        tmp3 = Input1[g]
        sum1 = tmp1.sum()
        sum3 = tmp3.sum()
        for i in tmp_pc:
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() 
            data3 = tmp3[start:end].sum()
            if data1 > 1.4 * data3:
                if g == '9' and i['start'] == 36600000 and i['end'] == 36800000:
                    print c
                pc[c].append((g , i['start'] , i['end']))
      
        
p = []      
for c in ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']:
    print c,len(pc[c]) / len(pc_data[c])   
    p.append(len(pc[c]) / len(pc_data[c]))
    
   
    

# Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
# fig = plt.figure(figsize = (14, 8))
# ax = fig.add_axes([Left  , HB , width , HH])
# ax.scatter(range(1,9) , p ,s = 50)
# ax.plot(range(1,9) , p)
# ax.set_xlim((0,9))
# ax.set_ylim((0,1))
# ax.set_xticks([1,2,3,4,5,6,7,8])
# ax.set_xticklabels(['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8'] ,fontsize = 20)
# ax.set_ylabel('Percentage' , fontsize = 30)
# run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_cluster_percentage_1.pdf')





Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (14, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.bar([1,2,3,4] , p[:4] , color = 'darkorange')
ax.bar([5,6,7,8] , p[4:] , color = 'steelblue')
ax.text(0.6 , p[0] + 0.05 , np.round(p[0] , 7))
ax.text(2 , p[1] + 0.05 , np.round(p[1] , 7))
ax.text(3 , p[2] + 0.05 , np.round(p[2] , 7))
ax.text(3.6 , p[3] + 0.05 , np.round(p[3] , 7))
ax.text(4.6 , p[4] + 0.05 , np.round(p[4] , 7))
ax.text(5.8 , p[5] + 0.05 , np.round(p[5] , 7))
ax.text(6.6 , p[6] + 0.05 , np.round(p[6] , 7))
ax.text(7.6 , p[7] + 0.05 , np.round(p[7] , 7))

ax.set_xlim((0,9))
ax.set_ylim((0,1))
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels(['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8'] ,fontsize = 20)
ax.set_ylabel('Percentage' , fontsize = 30)
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_cluster_percentage_barplot_1.pdf')





##------------------------H3K9me3_marked_compartment_related_genes------------------------------

gtf = Load_gtf('/public/home/lixinxin/data/ref/gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene_sites('/public/home/lixinxin/data/BDF1/RNA/RNA_new/mapping/allReps/FPKM/all_gene_expression.txt')



c1_gene =  Get_compartment_genes(pc['c1'] , union_gene)  
c4_gene =  Get_compartment_genes(pc['c4'] , union_gene)  
c5_gene =  Get_compartment_genes(pc['c5'] , union_gene)  
c6_gene =  Get_compartment_genes(pc['c6'] , union_gene)  
c7_gene =  Get_compartment_genes(pc['c7'] , union_gene)  
c8_gene =  Get_compartment_genes(pc['c8'] , union_gene)  


cl = {'A_To_B_Repro' : c1_gene , 'A_To_B_Resist' : c4_gene ,
      'B_To_A_Repro' : c5_gene , 'B_To_A_Partial' : c6_gene ,
      'B_To_A_Over' : c7_gene , 'B_To_A_Resist' : c8_gene }


head = ['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM']
for n in cl:
    OutputTXT(cl[n] , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_PC_related_genes/' + n + '_related_genes.txt' , head)
    
    
    gene_id = []
    for i in cl[n]:
        gene_name = i['gene_name']
        gene = gtf[gtf['gene_name'] == gene_name][0]
        gene_id.append(gene['gene_id'])
    
    out = open('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_PC_related_genes/' + n + '_related_genes_ID.txt' , 'w')    
    
    for i in gene_id:
        out.writelines(i + '\n')
    
    
    out.close()
    
    
   
    
##---------------------H3K9me3_marked_TAD_percent------------------------------------------

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
            start = i[1] // 100
            end = i[2] // 100
            data1 = tmp1[start:end].sum() 
            data3 = tmp3[start:end].sum()
            if data1 > 1.4 * data3:
                pc[c].append((g , i['start'] , i['end']))
      
        
p = []      
for c in ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']:
    print c,len(pc[c]) / len(pc_data[c])   
    p.append(len(pc[c]) / len(pc_data[c]))
    
   
    

Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (14, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.scatter(range(1,9) , p ,s = 50)
ax.plot(range(1,9) , p)
ax.set_xlim((0,9))
ax.set_ylim((0,1))
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels(['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8'] ,fontsize = 20)
ax.set_ylabel('Percentage' , fontsize = 30)
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_cluster_percentage.pdf')





Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (14, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.bar([1,2,3,4] , p[:4] , color = 'darkorange')
ax.bar([5,6,7,8] , p[4:] , color = 'steelblue')
ax.text(0.6 , p[0] + 0.05 , np.round(p[0] , 7))
ax.text(2 , p[1] + 0.05 , np.round(p[1] , 7))
ax.text(3 , p[2] + 0.05 , np.round(p[2] , 7))
ax.text(3.6 , p[3] + 0.05 , np.round(p[3] , 7))
ax.text(4.6 , p[4] + 0.05 , np.round(p[4] , 7))
ax.text(5.8 , p[5] + 0.05 , np.round(p[5] , 7))
ax.text(6.6 , p[6] + 0.05 , np.round(p[6] , 7))
ax.text(7.6 , p[7] + 0.05 , np.round(p[7] , 7))

ax.set_xlim((0,9))
ax.set_ylim((0,1))
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels(['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8'] ,fontsize = 20)
ax.set_ylabel('Percentage' , fontsize = 30)
run_Plot(fig , '/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/H3K9me3_marked_cluster_percentage_barplot.pdf')










CCS = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/CCS_Traditonal_PC_200K_Compartment_200K.txt' , skiprows = 0 , dtype = data_type )
n = 0
for g in chroms:
    tmp_pc = CCS[CCS['chr'] == g]
    tmp1 = Chip1[g]
    tmp3 = Input1[g]
    sum1 = tmp1.sum()
    sum3 = tmp3.sum()
    for i in range(len(tmp_pc)):
        start = (i * 200000) // 100
        end = (i + 1) * 200000 // 100
        data1 = tmp1[start:end].sum() 
        data3 = tmp3[start:end].sum()
        if data1 > data3:
            n += 1
            pc.append((g , start * 100 , end * 100))
  
    
  
    
  
    
