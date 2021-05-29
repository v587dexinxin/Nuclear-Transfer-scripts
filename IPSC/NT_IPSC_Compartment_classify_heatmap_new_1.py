# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:43:08 2021

@author: xxli
"""

from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from itertools import product

from matplotlib_venn import venn2, venn2_circles

# Our Own Color Map
# my_cmap = plt.get_cmap('bwr')
# my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['mediumblue' , 'white' , 'red'])
#my_cmap.set_bad('#2672a1')
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['U8' , np.float]})
NT_data_type = np.dtype({'names':['chr' , 'pos' , 'CCS' , 'NTs' , 'fESC'] , 
                         'formats':['U8' , np.int , np.float , np.float , np.float]})
IPSC_data_type = np.dtype({'names':['chr' , 'pos' , 'MEF' , 'IPSC_P3' , 'E14'] , 
                           'formats':['U8' , np.int , np.float , np.float , np.float]})
data_type = np.dtype({'names':['chr' , 'pos' ,'CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC_P3' , 'E14'] , 
                      'formats':['U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

lenth_type = np.dtype({'names':['chr' , 'len'] , 
                        'formats':['U8' , np.int]})
          
union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_type_2 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['U64' , 'U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
  


def get_union_gene_sites(Fil):
    '''
    '''
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F40_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if (i['CCS_FPKM'] > 0.1) or (i['NT5_FPKM'] > 0.1) or (i['NT6_FPKM'] > 0.1) or (i['F35_FPKM'] > 0.1) or (i['F40_FPKM'] > 0.1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS_FPKM'] ,  i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F40_FPKM']))
    union = np.array(union , dtype = union_type_1)
    return union


def Get_compartment_genes(pc , union_gene , Res):
    '''
    '''
    genes = []
    for i in pc:
        g = i[0]
        start = i[1] 
        end = i[1] +  Res
        union_tmp = union_gene[union_gene['chr'] == 'chr' + g]
        mask = (union_tmp['gene_site'] >= start) & (union_tmp['gene_site'] <= end)
        overlap = union_tmp[mask]
        if overlap.size != 0:
            for j in overlap:
                genes.append(j)
            
            
    genes = np.array(genes , dtype = union_gene.dtype)
    return genes



def get_raw_genes_new(Fil):
    '''
    '''
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT2' , 'NT3' , 'NT4' , 'NT5' , 'NT6' , 'F35' , 'F37' , 'F40' , 'F41'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 0.1) or (i['NT5'] > 0.1) or (i['NT6'] > 0.1) or (i['F35'] > 0.1) or (i['F40'] > 0.1):
            union.append(i)
    union = np.array(union , dtype = union_type)

    return union
    


def Load_gtf(gtfil):
    '''
    '''
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int , np.int]})
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


def Sort_c(c1 , c2):
    '''
    '''
    
    tmp1 = []
    tmp2 = []
    common = []
    for i in c1:
        if i not in c2:
            tmp1.append(i)
        else:
            common.append(i)
    for i in c2:
        if i not in c1:
            tmp2.append(i)
    c = common + tmp2
    n1 = sum([len(pc[x]) for  x in common])
    n2 = sum([len(pc[x]) for  x in tmp2])
    print ([keys.index(x) for x in c] , c , n1 , n2)
    return c , common , tmp2
            


def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()


def Sort(a , s_1):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1]))
    return a

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


CCS = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
CCS_1 = Sort(list(CCS) , 0)
CCS = np.array(CCS_1 , dtype = CCS.dtype)
NTs = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NTs_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
NTs = NTs[NTs['chr'] != 'X']
NTs_1 = Sort(list(NTs) , 0)
NTs = np.array(NTs_1 , dtype = CCS.dtype)
fESC = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\fESC_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']
fESC_1 = Sort(list(fESC) , 0)
fESC = np.array(fESC_1 , dtype = CCS.dtype)



MEF = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
MEF_1 = Sort(list(MEF) , 0)
MEF = np.array(MEF_1 , dtype = CCS.dtype)
IPS_P3 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
IPS_P3_1 = Sort(list(IPS_P3) , 0)
IPS_P3 = np.array(IPS_P3_1 , dtype = CCS.dtype)
IPS_P20 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P20_compartment_200K.txt' , dtype = pc_type)
IPS_P20 = IPS_P20[IPS_P20['chr'] != 'X']
IPS_P20_1 = Sort(list(IPS_P20) , 0)
IPS_P20 = np.array(IPS_P20_1 , dtype = CCS.dtype)
E14 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']
E14_1 = Sort(list(E14) , 0)
E14 = np.array(E14_1 , dtype = CCS.dtype)




union_pc = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_CCS)):
        union_pc.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc'] , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))

union_pc = np.array(union_pc , dtype = data_type)


##-----------------------------------classify-------------------------------------------------------------

c1 = ['A' , 'B']
c2 = ['A' , 'B']
c3 = ['A' , 'B']
c4 = ['A' , 'B']
c5 = ['A' , 'B']
c6 = ['A' , 'B']
loop_val = [c1,c2,c3,c4,c5,c6]
classify = []
for i in product(*loop_val):     
    classify.append(i[0] + i[1] + i[2] + i[3] + i[4] + i[5])


pc = {}

n = 0
for c in classify:
    if (c[0] == 'A' and c[1] == 'B' and c[2] == 'B' and c[-1] == 'B') or (c[0] == 'B' and c[1] == 'A' and c[2] == 'A' and c[-1] == 'A')  or (c[2] == 'B' and c[3] == 'A' and c[4] == 'B' and c[5] == 'B') or (c[2] == 'A' and c[3] == 'B' and c[4] == 'A' and c[5] == 'A'):
        pc[c] = []
    else:
        n += 1



for i in union_pc:
    k = ''
    data = [i['CCS'] , i['NTs'] , i['fESC'] , i['MEF'] , i['IPSC_P3'] , i['E14']]
    if data[0] > 0:
        k += 'A' 
    elif data[0] < 0:
        k += 'B'
    if data[1] > 0:
        k += 'A' 
    elif data[1] < 0:
        k += 'B'
    if data[2] > 0:
        k += 'A' 
    elif data[2] < 0:
        k += 'B'
    if data[3] > 0:
        k += 'A' 
    elif data[3] < 0:
        k += 'B'
    if data[4] > 0:
        k += 'A' 
    elif data[4] < 0:
        k += 'B'
    if data[5] > 0:
        k += 'A' 
    elif data[5] < 0:
        k += 'B'
    if k in pc.keys():
        pc[k].append(i)








##-----------------cluster_related_genes-------------------------------------------


mm10 = np.loadtxt('E:\\Data\\literature_data\\genome\\mm10_chr.txt' , dtype = lenth_type)
gtf = Load_gtf('E:\\Data\\literature_data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene_sites('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')


# IPSC_union_gene = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt' , dtype = union_type_1 , skiprows = 1)

# IPSC_union = []
# for i in IPSC_union_gene:
#     gene_site = int((i['start'] + i['end']) / 2)
#     IPSC_union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['F35'] , i['F40'] , i['MEF'] , i['IPS_P3'] , i['E14']))
# IPSC_union = np.array(IPSC_union , dtype = union_type_2)




genes = {}
for c in pc:
    if len(pc[c]) == 0:
        continue
    genes[c] = Get_compartment_genes(pc[c] , union_gene , 200000)
    print (c , len(pc[c]) , len(genes[c]))
    
    

    


    
    
    
##---------------------------------Specific_compartment_process-----------------------------------    
    
      
    
def plot_heatmap(matrix , vmin , vmax , pc , keys):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto' , interpolation = 'none')
    xticks = ['CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC' , 'E14']
    y = [0 ] + [len(pc[t]) for t in keys]
    yticks = [sum(y[:x]) for x in range(1 , len(y) + 1)]
    ax.set_xticks(np.arange(len(xticks)))
    ax.set_xticklabels(xticks,fontsize = 20)
    ax.set_yticks(yticks)
    ax.set_yticklabels([i + ':' + str(len(pc[i])) for i in  keys] + [''],fontsize = 0.5)
    # ax.set_ylabel('c1:' + str(len(pc['common_A_B'])) + ' , c2:' + str(len(pc['NT_A_B'])) + ' , c3:' + str(len(pc['IPSC_A_B'])) )
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Compartment numbers:' + str(len(matrix)),fontsize = 30)
    ##Colorbar
    ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
    cbar = fig.colorbar(im,cax = ax, orientation='vertical')
    cbar.set_ticks([vmin , vmax])
    
    return fig
    
 
    

k1 = ['ABBBBB' , 'ABBBAB' , 'ABBAAB']
k2 = ['BAAAAA' , 'BAAABA' , 'BAABBA']
k3 = ['BBBABB' , 'BABABB' , 'AABABB']
k4 = ['AAABAA' , 'ABABAA' , 'BBABAA']
k5 = ['ABBABB' , 'BAABAA']



    
classify = {'NT_A_B':k1 , 'NT_B_A':k2 , 'IPSC_A_B':k3 , 'IPSC_B_A':k4 , 'Common_repro':k5}

for cl in classify:
    n = 0 ; pc_new = {}
    for c in classify[cl]:
        pc_new[c] = pc[c]
        n += len(pc_new[c])
  
    matrix = np.zeros((n , 6))
    
    
    x = 0 
    for c in pc_new: 
        for i in range(x , x + len(pc_new[c])):
            matrix[i , 0] = pc_new[c][i - x]['CCS']
            matrix[i , 1] = pc_new[c][i - x]['NTs']
            matrix[i , 2] = pc_new[c][i - x]['fESC']
            matrix[i , 3] = pc_new[c][i - x]['MEF']
            matrix[i , 4] = pc_new[c][i - x]['IPSC_P3']
            matrix[i , 5] = pc_new[c][i - x]['E14']
        x += len(pc_new[c])
        print (c , len(pc_new[c]))

    
    # matrix_0 = Z_score(matrix)
    # plot_heatmap(matrix_0 , -2 ,1.5)
    fig = plot_heatmap(matrix , -0.06 , 0.06 , pc , classify[cl])
    run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\' + cl + '_compartment_heatmaplot.pdf')
    
            

    

##---------------------------------Specific_compartment_process_related_genes-----------------------------------           
        
classify_new = {'NT_A_B_Donor_specific':['ABBBBB' , 'ABBBAB'] , 'NT_A_B_Method_specific':['ABBAAB'],
                'IPSC_A_B_Donor_specific':['BBBABB' , 'BABABB'] , 'IPSC_A_B_Method_specific':['AABABB'],
                'NT_B_A_Donor_specific':['BAAAAA' , 'BAAABA'] , 'NT_B_A_Method_specific':['BAABBA'],
                'IPSC_B_A_Donor_specific':['AAABAA' , 'ABABAA'] , 'IPSC_B_A_Method_specific':['BBABAA']}
    

for C in classify_new:
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_Compartment_related_genes\\' + C + '_genes.txt' , 'w')
    out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_Compartment_related_genes\\' + C + '_genes_ID.txt' , 'w')
    out.writelines('\t'.join(['gene_name','chr','gene_site','CCS','NT5','NT6','F35','F40','MEF','IPS_P3','E14']) + '\n')
    for c in classify_new[C]:
        if c not in genes.keys():
            print (c)
            continue
        for i in genes[c]:
            out.writelines('\t'.join([str(x) for x in i]) + '\n')
            g = gtf[gtf['gene_name'] == i['gene_name']][0]
            out1.writelines(g['gene_id'] + '\n')
    out.close()
    out1.close()
    
        
    
    
    
    
##---------------------------------Chromatin_resistant_switch_common-----------------------------------       
 


c1 = ['A' , 'B']
c2 = ['A' , 'B']
c3 = ['A' , 'B']
c4 = ['A' , 'B']
c5 = ['A' , 'B']
c6 = ['A' , 'B']
loop_val = [c1,c2,c3,c4,c5,c6]
classify = []
for i in product(*loop_val):     
    classify.append(i[0] + i[1] + i[2] + i[3] + i[4] + i[5])


pc = {}

n = 0
for c in classify:
    if (c[0] == 'A' and c[1] == 'A' and c[2] == 'B' and c[-1] == 'B') or (c[0] == 'B' and c[1] == 'B' and c[2] == 'A' and c[-1] == 'A')  or (c[2] == 'B' and c[3] == 'A' and c[4] == 'A' and c[5] == 'B') or (c[2] == 'A' and c[3] == 'B' and c[4] == 'B' and c[5] == 'A'):
        pc[c] = []
    else:
        n += 1



for i in union_pc:
    k = ''
    data = [i['CCS'] , i['NTs'] , i['fESC'] , i['MEF'] , i['IPSC_P3'] , i['E14']]
    if data[0] > 0:
        k += 'A' 
    elif data[0] < 0:
        k += 'B'
    if data[1] > 0:
        k += 'A' 
    elif data[1] < 0:
        k += 'B'
    if data[2] > 0:
        k += 'A' 
    elif data[2] < 0:
        k += 'B'
    if data[3] > 0:
        k += 'A' 
    elif data[3] < 0:
        k += 'B'
    if data[4] > 0:
        k += 'A' 
    elif data[4] < 0:
        k += 'B'
    if data[5] > 0:
        k += 'A' 
    elif data[5] < 0:
        k += 'B'
    if k in pc.keys():
        pc[k].append(i)









k1 = ['AABBBB' , 'AABBAB' , 'AABABB']
k2 = ['BBAAAA' , 'BBAABA' , 'BBABAA']
k3 = ['BBBAAB' , 'BABAAB' , 'ABBAAB']
k4 = ['AAABBA' , 'ABABBA' , 'BAABBA']
k5 = ['AABAAB'] 
k6 = ['BBABBA']


    
classify = {'Resist' : k1 + k3 + k5 +  k2 + k4 + k6 }

for cl in classify:
    n = 0 ; pc_new = {}
    for c in classify[cl]:
        pc_new[c] = pc[c]
        n += len(pc_new[c])
  
    matrix = np.zeros((n , 6))
    
    
    x = 0 
    for c in pc_new: 
        for i in range(x , x + len(pc_new[c])):
            matrix[i , 0] = pc_new[c][i - x]['CCS']
            matrix[i , 1] = pc_new[c][i - x]['NTs']
            matrix[i , 2] = pc_new[c][i - x]['fESC']
            matrix[i , 3] = pc_new[c][i - x]['MEF']
            matrix[i , 4] = pc_new[c][i - x]['IPSC_P3']
            matrix[i , 5] = pc_new[c][i - x]['E14']
        x += len(pc_new[c])
        print (c , len(pc_new[c]))

    
    # matrix_0 = Z_score(matrix)
    # plot_heatmap(matrix_0 , -2 ,1.5)
    fig = plot_heatmap(matrix , -0.06 , 0.06 , pc , classify[cl])
    run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\Resistant_compartment\\' + cl + '_compartment_heatmaplot.pdf')
    
            

##---------------------------------resistant_compartment_process_related_genes-----------------------------------    


       
genes = {}
for c in pc:
    if len(pc[c]) == 0:
        continue
    genes[c] = Get_compartment_genes(pc[c] , union_gene , 200000)
    print (c , len(pc[c]) , len(genes[c]))
    
        
classify_new = {'NT_A_B_resistant':k1 , 'NT_B_A_resistant':k2,
                'IPSC_A_B_resistant':k3 , 'IPSC_B_A_resistant':k4,
                'NT_IPSC_A_B_common':k5 , 'NT_IPSC_B_A_common':k6}
    

for C in classify_new:
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_resistant_Compartment_related_genes\\' + C + '_genes.txt' , 'w')
    out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\NT_IPSC_compartment_classify\\IPSC_P3\\NT_IPSC_P3_resistant_Compartment_related_genes\\' + C + '_genes_ID.txt' , 'w')
    out.writelines('\t'.join(['gene_name','chr','gene_site','CCS','NT5','NT6','F35','F40','MEF','IPS_P3','E14']) + '\n')
    for c in classify_new[C]:
        if c not in genes.keys():
            print (c)
            continue
        for i in genes[c]:
            out.writelines('\t'.join([str(x) for x in i]) + '\n')
            g = gtf[gtf['gene_name'] == i['gene_name']][0]
            out1.writelines(g['gene_id'] + '\n')
    out.close()
    out1.close()
    



# n = 0 ; pc_new = {} ; keys_new = []
# for c in pc:
#     c_new = c[0] + c[2] + c[4] + c[1] + c[3] + c[5]
#     pc_new[c_new] = np.array(pc[c] ,  dtype=data_type)
#     n += len(pc_new[c_new])
#     keys_new.append(c_new)

     

# keys = ['AABAAA' , 'AABABA' , 'AABABB' , 'AABBBB' , 'AABBBA' , 'AABBAA' , 'AABBAB' , 'AAAAAB' , 'ABAAAB' , 'ABBAAB' , 'BBBAAB' , 'BBAAAB' , 'BAAAAB' , 'BABAAB' , 'AABAAB' , 
#           'BBAAAA' , 'BBAAAB' , 'BBAABA' , 'BBAABB' , 'BBABBB' , 'BBABAA' , 'BBABAB' , 'AAABBA' , 'AABBBA' , 'ABABBA' , 'ABBBBA' , 'BBBBBA' , 'BAABBA' , 'BABBBA' , 'BBABBA']
    
# pc_resis = {};tmp = [] ; n = 0 ; keys_new = []
# for c in keys:
#     if len(pc_new[c]) == 0:
#         print (c)
#     pc_resis[c] = pc_new[c]
#     tmp.extend(list(pc_new[c]))
#     n += len(pc_resis[c])
#     keys_new.append(c)
# tmp = np.array(tmp , dtype = pc_resis['AABABB'].dtype)



# matrix = np.zeros((len(tmp) , 6))

# for i in range(len(tmp)):
#     matrix[i , 0] = tmp[i]['CCS']
#     matrix[i , 1] = tmp[i]['NTs']
#     matrix[i , 2] = tmp[i]['fESC']
#     matrix[i , 3] = tmp[i]['MEF']
#     matrix[i , 4] = tmp[i]['IPSC_P3']
#     matrix[i , 5] = tmp[i]['E14']
    

    
    
    
# fig = plot_heatmap(matrix , -0.08 , 0.06 , pc_resis , keys_new)    
# run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S7\\NT_IPSC_chromatin_resistant_compartment.pdf')   
        


# ##-------------------------------------A-B-----------------------------------------------
# Ab = 0 ; aB = 0 ; AB = 0

# NT_A_B = np.hstack([pc_resis[c] for c in ['AABAAA' , 'AABABA' , 'AABABB' , 'AABBBB' , 'AABBBA' , 'AABBAA' , 'AABBAB']])
# IPSC_A_B = np.hstack([pc_resis[c] for c in ['AAAAAB' , 'ABAAAB' , 'ABBAAB' , 'BBBAAB' , 'BBAAAB' , 'BAAAAB' , 'BABAAB']])
# Common_A_B = pc_resis['AABAAB']
# NT_B_A = np.hstack([pc_resis[c] for c in ['BBAAAA' , 'BBAAAB' , 'BBAABA' , 'BBAABB' , 'BBABBB' , 'BBABAA' , 'BBABAB']])
# IPSC_B_A = np.hstack([pc_resis[c] for c in ['AAABBA' , 'AABBBA' , 'ABABBA' , 'ABBBBA' , 'BBBBBA' , 'BAABBA' , 'BABBBA']])
# Common_B_A = pc_resis['BBABBA']




# Ab = len(NT_A_B)  
# aB = len(IPSC_A_B)   
# AB = len(Common_A_B)
    


# print (Ab , aB , AB)

# fig = plt.figure(figsize = (10, 10))
# venn2(subsets = (Ab , aB , AB), set_labels = ('NT_process', 'IPSC_process'))
# run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S7\\NT_IPSC_chromatin_resistant_compartment_A_B_Venn2.pdf')


# ##-------------------------------------B-A-----------------------------------------------

 


# Ab = len(NT_B_A)  
# aB = len(IPSC_B_A)   
# AB = len(Common_B_A)
    

# print (Ab , aB , AB)

# fig = plt.figure(figsize = (10, 10))
# venn2(subsets = (Ab , aB , AB), set_labels = ('NT_process', 'IPSC_process'))
# run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S7\\NT_IPSC_chromatin_resistant_compartment_B_A_Venn2.pdf')






# compaerment = {'NT_speci_A_TO_B_resis':NT_A_B , 'IPSC_speci_A_TO_B_resis':IPSC_A_B , 'NT_IPSC_A_TO_B_common_resis':Common_A_B , 
#                'NT_speci_B_TO_A_resis':NT_B_A , 'IPSC_speci_B_TO_A_resis':IPSC_B_A , 'NT_IPSC_B_TO_A_common_resis':Common_B_A}    
    

# genes = {} ; gene_id = {}

# for c in compaerment:
#     genes[c] = Get_compartment_genes(compaerment[c] , union_gene , 200000)
    
    
# for c in genes:
#     out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\compartment_resist_gene\\' + c + 'compartment_related_genes.txt' , 'w')
#     out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\compartment_resist_gene\\' + c + 'compartment_related_genes_ID.txt' , 'w')
#     out.writelines('\t'.join(list(genes[c].dtype.names)) + '\n')
#     for i in genes[c]:
#         out.writelines('\t'.join([str(x) for x in i]) + '\n')
#         g = gtf[gtf['gene_name'] == i['gene_name']][0]
#         out1.writelines(g['gene_id'] + '\n')
#     out.close()
#     out1.close()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    