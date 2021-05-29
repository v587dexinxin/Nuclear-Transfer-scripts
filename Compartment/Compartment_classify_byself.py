# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 21:00:14 2020

@author: xxli
"""


from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
import pandas as pd
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
# from palettable.colorbrewer.qualitative import  Paired_10
# color = Paired_10.hex_colors

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['mediumblue' , 'white' , 'red'])
#my_cmap.set_bad('#2672a1')
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['U8' , np.float]})
                    
 

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

                   
def Classify3_compartment(comp1 , comp2 , comp3):
    classify = {'AAA':[] , 'AAB':[] , 'ABA':[] , 'ABB':[] , 'BBB':[] , 'BBA':[] , 'BAA':[] , 'BAB':[]}
    for g in chrom:
        tmp_comp1 = comp1[comp1['chr'] == g]
        tmp_comp2 = comp2[comp2['chr'] == g]
        tmp_comp3 = comp3[comp3['chr'] == g]
        for i in range(len(tmp_comp1)):
            score1 = tmp_comp1[i]['pc']
            score2 = tmp_comp2[i]['pc']
            score3 = tmp_comp3[i]['pc']
            if (score1 > 0) and (score2 > 0) and (score3 > 0):
                classify['AAA'].append((score1 , score2 , score3 , g , i))
            elif (score1 > 0) and (score2 > 0) and (score3 < 0):
                classify['AAB'].append((score1 , score2 , score3 , g , i))
            elif (score1 > 0) and (score2 < 0) and (score3 > 0):
                classify['ABA'].append((score1 , score2 , score3 , g , i))
            elif (score1 > 0) and (score2 < 0) and (score3 < 0):
                classify['ABB'].append((score1 , score2 , score3 ,g , i))
            elif (score1 < 0) and (score2 < 0) and (score3 < 0):
                classify['BBB'].append((score1 , score2 , score3 , g , i))
            elif (score1 < 0) and (score2 < 0) and (score3 > 0):
                classify['BBA'].append((score1 , score2 , score3 , g , i))
            elif (score1 < 0) and (score2 > 0) and (score3 > 0):
                classify['BAA'].append((score1 , score2 , score3 , g , i))
            elif (score1 < 0) and (score2 > 0) and (score3 < 0):
                classify['BAB'].append((score1 , score2 , score3 , g , i))
    return classify


def get_3cell_matrix(matrix_raw):
    matrix = np.zeros((len(matrix_raw) , 3))
    pos = []
    n = 0
    pos.append(('chr' , 'start' , 'end' , 'CCS_pc1' , 'NTs_pc1' , 'fESC_pc1'))
    for i in matrix_raw:
        g = i[-2]
        start = i[-1] * 200000
        end = (i[-1] + 1) * 200000
        matrix[n , 0] = i[0]
        matrix[n , 1] = i[1]
        matrix[n , 2] = i[2]
        pos.append((g , start , end , i[0] , i[1] , i[2]))
        n += 1
    return matrix , pos
    



def plot_heatmap(matrix , vmin , vmax):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    x = ['CCS','ntS','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 10)
    ax.set_ylabel('c1:' + str(len(Matrix1)) + ' , c2:' + str(len(Matrix2)) + ' , c3:' + str(len(Matrix3)) + ' , c4:' + str(len(Matrix4)) + ' , c5:' + str(len(Matrix5)) + ' , c6:' + str(len(Matrix6)) + ' , c7:' + str(len(Matrix7)) + ' , c8:' + str(len(Matrix8)))
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Compartment numbers:' + str(len(matrix)),fontsize = 30)
    ##Colorbar
    ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
    cbar = fig.colorbar(im,cax = ax, orientation='vertical')
    cbar.set_ticks([vmin , vmax])
    
    return fig



def get_union_gene_sites(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
         
                 
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
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene
    


def Load_gtf(gtfil):
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



        
    
CCS = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT5_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT6_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
F35 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F35_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
F35 = F35[F35['chr'] != 'X']
F40 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F40_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
F40 = F40[F40['chr'] != 'X']


NTs = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NTs_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
NTs = NTs[NTs['chr'] != 'X']
fESC = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\fESC_Traditonal_PC_200K_Compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']

gtf = Load_gtf('E:\\Data\\literature_data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
union_gene = get_union_gene_sites('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')
union_gene_raw = get_raw_genes_new('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')



Res = 200000
chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']


classify = Classify3_compartment(CCS , NTs , fESC)


matrix1 = [] ; matrix2 = [] ; matrix3 = [] ; matrix4 = [] ; matrix5 = [] ; matrix6 = [] ; matrix7 = [] ; matrix8 = []

Sort(classify['ABB'] , 0 , 0) ; Sort(classify['BAA'] , 0 , 0)
for i in classify['ABB']:
    if i[1] > 0.2 * i[2]:
        matrix1.append(i)
    elif i[1] < 1.5 * i[2]:
        matrix2.append(i)
    else:
        matrix3.append(i)
        
matrix4 = Sort(classify['AAB'] , 0 , 0)        

for i in classify['BAA']:
    if i[1] < 0.2 * i[2]:
        matrix5.append(i)
        
    elif i[1] > 1.5 * i[2]:
        matrix6.append(i)
    else:
        matrix7.append(i)
        
matrix8 = Sort(classify['BBA'] , 0 , 0)
        
Matrix1,pos1 = get_3cell_matrix(matrix3) ; Matrix2,pos2 = get_3cell_matrix(matrix1)
Matrix3,pos3 = get_3cell_matrix(matrix2) ; Matrix4,pos4 = get_3cell_matrix(matrix4) 
Matrix5,pos5 = get_3cell_matrix(matrix7) ; Matrix6,pos6 = get_3cell_matrix(matrix5)
Matrix7,pos7 = get_3cell_matrix(matrix6) ; Matrix8,pos8 = get_3cell_matrix(matrix8)
 

          
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster1_pc1.txt' , pos1)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster2_pc1.txt' , pos2)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster3_pc1.txt' , pos3)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster4_pc1.txt' , pos4)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster5_pc1.txt' , pos5)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster6_pc1.txt' , pos6)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster7_pc1.txt' , pos7)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\Compartment_cluster8_pc1.txt' , pos8)


#-------------------------plot_compartment_states-------------------------------
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['mediumblue' , 'orange'])
my_cmap.set_bad('#2672a1')

matrix = np.array([[1 , -1 , -1] , [1 , -1 , -1] , [1 , -1 , -1] , [1 , 1 , -1] , [-1 , 1 , 1] , [-1 , 1 , 1] , [-1 , 1 , 1] , [-1 , -1 , 1]])
left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (6, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
im = ax.imshow(matrix,vmax=1,vmin=-1,cmap=my_cmap,aspect = 'auto')
x = ['CCS','NT','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 10)
ax.set_yticks(np.arange(len(matrix)))
ax.set_yticklabels(['ABrB' , 'ABpB' , 'ABoB' , 'AAB' , 'BArA' , 'BApA' , 'BAoA' , 'BBA'],fontsize = 10)
ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
for i in range(len(matrix)):
    ax.plot([0 - 0.5 , 3 - 0.5] , [i - 0.5 , i - 0.5] , c = 'black')
for j in range(len(x)):
    ax.plot([j - 0.5 , j - 0.5] , [0 - 0.5 , 8 - 0.5] , c = 'black')
        
plt.title('Compartment states',fontsize = 30)

    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_states.pdf')
      
 
#-------------------------plot_compartment_states_number-------------------------------               
c1 = len(Matrix1) / 546 ; c2 = len(Matrix2) / 546 ; c3 = len(Matrix3) / 546 ; c4 = len(Matrix4) / 546
c5 = len(Matrix5) / 1889 ; c6 = len(Matrix6) / 1889 ; c7 = len(Matrix7) / 1889 ; c8 = len(Matrix8) / 1889

left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (6, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))

ax.bar(1 , c1 + c2 + c3 + c4 , color='green')
ax.bar(1 , c1 + c2 + c3 , color='red')
ax.bar(1 , c1 + c2 , color='orange')
ax.bar(1 , c1 , color='mediumblue')

ax.bar(2 , c5 + c6 + c7 + c8 , color='green')
ax.bar(2 , c5 + c6 + c7 , color='red')
ax.bar(2 , c5 + c6 , color='orange')
ax.bar(2 , c5 , color='mediumblue')



x = ['A-B','B-A']
ax.set_xticks([1 , 2])
ax.set_xticklabels(x,fontsize = 10)
ax.set_ylabel('Percent(%)',fontsize = 10)
ax.set_ylim((0.5 , 2.5))
ax.set_ylim((0 , 1))

        
plt.title('Compartment states numbers',fontsize = 10)

    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_states_numbers.pdf') 




#--------------------Compartment_classify-----------------------------------

cl = {'A_To_B_Reprogramming_compartment' : pos1 , 'A_To_B_Partial_compartment' : pos2 , 
      'A_To_B_Over_compartment' : pos3 , 'A_To_B_Resistant_compartment' : pos4 ,
      'B_To_A_Reprogramming_compartment' : pos5 , 'B_To_A_Partial_compartment' : pos6 ,
      'B_To_A_Over_compartment' : pos7 , 'B_To_A_Resistant_compartment' : pos8 }

for n in cl:
    out_fil = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\fc_0.2_1.5\\' + n + '_PC1_nofiltering.txt'
    Write2fils_nochr(out_fil , cl[n])
    
    
    


#--------------------Compartment_genes--------------------------------------

c1_gene =  Get_compartment_genes(pos1[1:] , union_gene , Res)  
c2_gene =  Get_compartment_genes(pos2[1:] , union_gene , Res)  
c3_gene =  Get_compartment_genes(pos3[1:] , union_gene , Res)  
c4_gene =  Get_compartment_genes(pos4[1:] , union_gene , Res)  
c5_gene =  Get_compartment_genes(pos5[1:] , union_gene , Res)  
c6_gene =  Get_compartment_genes(pos6[1:] , union_gene , Res)  
c7_gene =  Get_compartment_genes(pos7[1:] , union_gene , Res)  
c8_gene =  Get_compartment_genes(pos8[1:] , union_gene , Res)  


cl = {'A_To_B_Reprogramming_compartment' : c1_gene , 'A_To_B_Partial_compartment' : c2_gene , 
      'A_To_B_Over_compartment' : c3_gene , 'A_To_B_Resistant_compartment' : c4_gene ,
      'B_To_A_Reprogramming_compartment' : c5_gene , 'B_To_A_Partial_compartment' : c6_gene ,
      'B_To_A_Over_compartment' : c7_gene , 'B_To_A_Resistant_compartment' : c8_gene }


head = ['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM']
for n in cl:
    OutputTXT(cl[n] , 'F:\\work\\ntESC_3Dreprogramming\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\' + n + '_related_genes_nofiltering.txt' , head)
    
    
    gene_id = []
    for i in cl[n]:
        gene_name = i['gene_name']
        gene = gtf[gtf['gene_name'] == gene_name][0]
        gene_id.append(gene['gene_id'])
    
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\' + n + '_related_genes_ID_nofiltering.txt' , 'w')    
    
    for i in gene_id:
        out.writelines(i + '\n')
    
    
    out.close()


head = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCs_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F40_FPKM']
for n in cl:
        
    Gene = []
    for i in cl[n]:
        gene_name = i['gene_name']
        gene = gtf[gtf['gene_name'] == gene_name][0]
        gene1 = union_gene_raw[union_gene_raw['gene_name'] == gene_name][0]
        Gene.append([gene['gene_id'] , i['gene_name'] , i['chr'] , gene1['strand'] , gene1['start'] , gene1['end'] , i['CCS_FPKM'] , i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F40_FPKM']])
    OutputTXT(Gene , 'F:\\work\\ntESC_3Dreprogramming\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\excel\\' + n + '_related_genes_nofiltering.txt' , head)


    
    
#--------------------compartment_related_Diff_genes-------------------------------------- 
    
diff_gene = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\diff_expression_q_0.05\\RNA_fESC_vs_NTs_diff_expression_q_0.05.csv' , index_col=0)
gene_names = []
for i in diff_gene.index:
    g_name = gtf[gtf['gene_id'] == i][0]['gene_name']
    gene_names.append(g_name)
    
data_type = np.dtype({'names':['gene_name' , 'diff'],
                      'formats':['U64' , np.float]})

diff = {'A_To_B_Reprogramming_compartment' : [] , 'A_To_B_Partial_compartment' : [] , 
      'A_To_B_Over_compartment' : [] , 'A_To_B_Resistant_compartment' : [] ,
      'B_To_A_Reprogramming_compartment' : [] , 'B_To_A_Partial_compartment' : [] ,
      'B_To_A_Over_compartment' : [] , 'B_To_A_Resistant_compartment' : []}

for c in cl:
    for i in cl[c]:
        if i['gene_name'] in gene_names:
            diff[c].append(i)

head = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCs_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F40_FPKM' , 'log2FoldChange']
for n in diff:     
    Gene = []
    print (n , len(diff[n]))
    for i in diff[n]:
        gene_name = i['gene_name']
        gene = gtf[gtf['gene_name'] == gene_name][0]
        gene1 = union_gene_raw[union_gene_raw['gene_name'] == gene_name][0]
        Gene.append([gene['gene_id'] , i['gene_name'] , i['chr'] , gene1['strand'] , gene1['start'] , gene1['end'] , i['CCS_FPKM'] , i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F40_FPKM'] , diff_gene.loc[gene['gene_id']]['log2FoldChange']])
    OutputTXT(Gene , 'F:\\work\\ntESC_3Dreprogramming\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\Compartment_classify_byself\\genes\\fc_0.2_1.5\\Compartment_related_DEGs\\DEG_q_0.05\\' + n + '_related_genes_nofiltering.txt' , head)

            
            
            
            
            

#--------------------Diff_genes--------------------------------------
    
diff_gene = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\diff_expression_q_0.05_fc_1.5\\ESC_vs_NTs_diff_expression_q_0.05_fc1.5.csv' , index_col=0)

    
data_type = np.dtype({'names':['gene_name' , 'diff'],
                      'formats':['U64' , np.float]})

diff = {'A_To_B_Reprogramming_compartment' : [] , 'A_To_B_Partial_compartment' : [] , 
      'A_To_B_Over_compartment' : [] , 'A_To_B_Resistant_compartment' : [] ,
      'B_To_A_Reprogramming_compartment' : [] , 'B_To_A_Partial_compartment' : [] ,
      'B_To_A_Over_compartment' : [] , 'B_To_A_Resistant_compartment' : []}

for c in cl:
    for i in cl[c]:
        if i['gene_name'] in diff_gene.index:
             if (i['NT5_FPKM'] > i['F35_FPKM']) and (i['NT5_FPKM'] > i['F40_FPKM']) and (i['NT6_FPKM'] > i['F35_FPKM']) and (i['NT6_FPKM'] > i['F40_FPKM']):
                nts = (i['NT5_FPKM'] + i['NT6_FPKM']) / 2
                fesc = (i['F35_FPKM'] + i['F40_FPKM']) / 2
                d = nts - fesc
                diff[c].append((i['gene_name'] , fesc))
        
    diff[c].sort(key = lambda x:x[1])
    diff[c] = np.array(diff[c] , dtype = data_type)
    
    
        
        
gene_list = list(diff['A_To_B_Partial_compartment']['gene_name'][-5:])
gene_list = list(diff['A_To_B_Resistant_compartment']['gene_name'][-5:])
gene_list = list(diff['B_To_A_Over_compartment']['gene_name'][-5:])
gene_list.reverse()


fig = plt.figure(figsize = (12, 12))
for i , gene in enumerate(gene_list):
    tmp = union_gene[union_gene['gene_name'] == gene][0]
    tmp_ccs = tmp['CCS']
    tmp_nt2 = tmp['NT2']
    tmp_nt3 = tmp['NT3']
    tmp_nt4 = tmp['NT4']
    tmp_nt5 = tmp['NT5']
    tmp_nt6 = tmp['NT6']
    tmp_f35 = tmp['F35']
    tmp_f37 = tmp['F37']
    tmp_f40 = tmp['F40']
    tmp_f41 = tmp['F41']
    
    y = [tmp_ccs , tmp_nt2 , tmp_nt3 , tmp_nt4 , tmp_nt5 , tmp_nt6 , tmp_f35 , tmp_f37 , tmp_f40 , tmp_f41]
    plt.plot(range(len(y)),[np.log2(x + 1) for x in y], c = color[i], marker = 'o', label = gene)
    plt.xticks(range(10), ['CCS','NT2','NT3','NT4','NT5','NT6','F35','F37','F40','F41'])
    plt.legend(loc = 'upper right')
    plt.ylabel('Log2(fpkm+1)')
    plt.title('B_to_A_Over')


run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\Plot\B_to_A_Over_max_diff_gene_FPKM_linechart_all.pdf')


    
    
    
    
fig = plt.figure(figsize = (12, 12))
for i , gene in enumerate(gene_list):
    tmp = union_gene[union_gene['gene_name'] == gene][0]
    tmp_ccs = tmp['CCS']
    tmp_nt2 = tmp['NT2']
    tmp_nt3 = tmp['NT3']
    tmp_nt4 = tmp['NT4']
    tmp_nt5 = tmp['NT5']
    tmp_nt6 = tmp['NT6']
    tmp_f35 = tmp['F35']
    tmp_f37 = tmp['F37']
    tmp_f40 = tmp['F40']
    tmp_f41 = tmp['F41']
    
    y = [tmp_ccs , tmp_nt5 , tmp_nt6 , tmp_f35 , tmp_f40]
    plt.plot(range(len(y)),[np.log2(x + 1) for x in y], c = color[i], marker = 'o', label = gene)
    plt.xticks(range(len(y)), ['CCS','NT5','NT6','F35','F40'])
    plt.legend(loc = 'upper right')
    plt.ylabel('Log2(fpkm+1)')
    plt.title('B_to_A_Over')


run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\Plot\\B_to_A_Over_max_diff_gene_FPKM_linechart.pdf')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    