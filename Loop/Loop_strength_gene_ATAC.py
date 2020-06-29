# -*- coding: utf-8 -*-
"""
Created on Tue Nov 06 16:34:10 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import sys, cPickle
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
from sklearn.cluster import KMeans
from scipy import stats

## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#D3D3D3')
def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a


def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)
    
    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1


def Load_loop_strength(LoopFil):
    union_type = np.dtype({'names':['chr','start','end','CCS','NT5','NT6','fESC'],
                       'formats':['S2' , np.int , np.int , np.float , np.float , np.float , np.float]})
    
    union_sites = np.loadtxt(LoopFil , skiprows = 1 , dtype = union_type)        
    sort_union_sites = sorted(union_sites , key = lambda x:(x[0],x[1]))
    union_sites = np.array(sort_union_sites , dtype = union_type)
    return union_sites


def Load_loop_strength_in_pairs(union_sites):
    al = {}
    ccs = []
    nt5 = []
    nt6 = []
    fesc = []
    loop = []
    for i in union_sites:
        loop.append((i['chr'] , i['start'] , i['end']))
        ccs.append(i['CCS'])
        nt5.append(i['NT5'])
        nt6.append(i['NT6'])
        fesc.append(i['fESC'])
    
    
    for i in range(len(ccs) - 1 , -1 , -1):
        if (ccs[i] > 50) or (nt5[i] > 50) or (nt6[i] > 50) or (fesc[i] > 50):

#        print i
            del ccs[i]
            del nt5[i]
            del nt6[i]
            del fesc[i]
            del loop[i]
    
    al = {'CCS_NT5':[ccs,nt5] , 'CCS_NT6':[ccs,nt6] , 'CCS_fESC':[ccs,fesc] , 'NT5_NT6':[nt5,nt6] , 'NT5_fESC':[nt5,fesc] , 'NT6_fESC':[nt6,fesc]}
    return ccs , nt5 , nt6 , fesc , loop , al

def union_loop(a,peak_type):
    '''
    a: np.array of peaks ,dtype = peak_type
    '''
    peak_new = []
    a = np.array(a , dtype = peak_type)
    for i in a:
        s = i['start']
        e = i['end']
        mask = (s == a[a['chr'] == i['chr']]['start']) & (e == a[a['chr'] == i['chr']]['end'])
        overlap = a[a['chr'] == i['chr']][mask]
        if overlap.size != 0:
            peak_new.append(i)
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type)
        else:
            continue
    
    peak_new = Sort(peak_new , 0 ,1)
    peak_new = np.array(peak_new , dtype = peak_type)
    return peak_new



def Select_diff_loop(loop , al):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    vs_loop = {}    
    for c in al:
        vs_loop[c] = []
        for i in range(len(al[c][0])):
            l1 = al[c][0][i]
            l2 = al[c][1][i]
            vs_loop[c].append((loop[i][0] , float(loop[i][1]) , float(loop[i][2]) , l1 , l2 ))
    
    diff_loop = []
    for c in vs_loop:
        for i in vs_loop[c]:
            if (i[-1] >= 3) | (i[-2] >= 3) and ((i[-1] >= 2 * i[-2]) or (i[-1] <= 0.5 * i[-2])):
                diff_loop.append((i[0] , float(i[1]) , float(i[2])))
    
    diff_loop = union_loop(diff_loop , loop_type)   
    return diff_loop



def Heat_map(matrix , title):            
    left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    #matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    #matrix = matrix_1[:,:-1]
    vmax = np.percentile(matrix,90)
    vmin = np.percentile(matrix,10)
    im = ax.imshow(matrix , vmax = vmax , vmin = vmin , cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 15)
    plt.title(title,fontsize = 20)


def center_sites(lists):
    sum_x = 0
    sum_y = 0
    for x in lists:
        sum_x += x[1]
        sum_y += x[2]
    n = len(lists)
    return [float(sum_x)/n, float(sum_y)/n]



def union_peaks(peaks):
    peak_type = ({'names':['chr' , 'start' , 'end' , 'score' , 'cell'],
                  'formats':['S8' , np.int , np.int , np.float , 'S8']})    
    chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
    union_type = ({'names':['chr' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                  'formats':['S8' , np.int , np.int , np.float , np.float , np.float , np.float]})    
    
    peaks_new = []
    for c in cell:
        for i in peaks[c]:
            peaks_new.append((i['chr'] , i['start'] , i['end'] , i['score'] , c))
    sort_peaks = Sort(peaks_new , 0 , 1)
    sort_peaks = np.array(sort_peaks , dtype = peak_type)
    union = []
    tmp = []
    for g in chrom:
        lib = sort_peaks[sort_peaks['chr'] == 'chr' + g]
        lib1 = sort_peaks[sort_peaks['chr'] == 'chr' + g]
        for i in lib:
            start = i['start']
            end = i['end']
            mask = (start <= lib1['end']) & (end >= lib1['start'])
            overlap = lib1[mask]
            if overlap.size > 1:
                if i in tmp:
                    pass
                else:
                    if tmp:
                        center = center_sites(tmp)
                        tmp1 = {'CCS' : 0 , 'NT5' : 0 , 'NT6' : 0 , 'fESC' : 0}
                        for k in tmp:
                            tmp1[k['cell']] = k['score']
                        union.append((g  , int(center[0]) , int(center[1]) , tmp1['CCS'] , tmp1['NT5'] , tmp1['NT6'] , tmp1['fESC']))
                    else:
                        pass
                    
                    tmp = []
                lib1 = list(lib1)
                for j in overlap:
                    tmp.append(j)
                    lib1.remove(j)
                lib1 = np.array(lib1 , dtype = peak_type)
            
            else:
                if i in tmp:
                    pass
                else:
                    if tmp:
                        center = center_sites(tmp)
                        tmp1 = {'CCS' : 0 , 'NT5' : 0 , 'NT6' : 0 , 'fESC' : 0}
                        for k in tmp:
                            tmp1[k['cell']] = k['score']
                        union.append((g  , int(center[0]) , int(center[1]) , tmp1['CCS'] , tmp1['NT5'] , tmp1['NT6'] , tmp1['fESC']))
                    else:
                        pass
                    
                    tmp = []
                    tmp.append(i)
                    lib1 = list(lib1)
                    lib1.remove(i)
                    lib1 =  np.array(lib1 , dtype = peak_type)
    
    union = np.array(union , dtype = union_type)
    return union


LoopFil = 'D:\\workspace\\loop\\loop_strength_20K\\union_sites_ave_20K.txt'
union_sites = Load_loop_strength(LoopFil)
(ccs , nt5 , nt6 , fesc , loop , al) = Load_loop_strength_in_pairs(union_sites)


#-------------------------------------------matrix-------------------------------------------------------------------------

diff_loop = Select_diff_loop(loop , al)


matrix = np.zeros((len(diff_loop) , 5))
for i in range(len(matrix)):
    index = loop.index((diff_loop[i][0] , diff_loop[i][1] , diff_loop[i][2]))
    matrix[i][0] = np.log2(ccs[index])
    matrix[i][1] = np.log2(nt5[index])
    matrix[i][2] = np.log2(nt6[index])
    matrix[i][3] = np.log2(fesc[index])
    matrix[i][4] = index


##--------------------------------------------Z-score 标准化------------------------------------------------------------------------   
matrix_1 = stats.zscore(matrix[: , :-1] , axis = 1 , ddof = 1)
matrix_1[np.isnan(matrix_1)] = 0
matrix = np.hstack((matrix_1 , matrix[:,4].reshape((len(matrix) , 1))))

matrix_1 = K_means_cluster(matrix , matrix[:,:-1] , 8)  

#------------------------------------------Heatmap_Plot----------------------------------------------------------------------------

Heat_map(matrix_1[:,:-2] , 'diff_loop_strength loop numbers:' + str(len(matrix)))


#----------------------------------------adjust_classify---------------------------------------------------------------------------
loop_type = ({'names':['chr' , 'start' , 'end'],
              'formats':['S8' , np.int , np.int]})
classify_matrix = {}
for i in matrix_1:
    keys = i[-1]
    if keys not in classify_matrix.keys():
        classify_matrix[keys] = i
    else:
        classify_matrix[keys] = np.vstack((classify_matrix[keys] , i))


matrix_new =  classify_matrix[1]   
for i in [3 , 7 , 6 , 0 , 5 , 2 , 4]:
    matrix_new = np.vstack((matrix_new , classify_matrix[i]))


Heat_map(matrix_new[:,:-2] , 'diff_loop_strength loop numbers:' + str(len(matrix)))    

index = matrix_new[: , -2:]
index
loop_new = {}
for i in index:
    if i[1] not in loop_new.keys():
        loop_new[i[1]] = []
        loop_new[i[1]].append(loop[int(i[0])])
    else:
        loop_new[i[1]].append(loop[int(i[0])])


loop_new
loop_new = {}
for i in range(8):
    loop_new[i] = []

for i in index:
    if i[1] == 1:
        loop_new[0].append(loop[int(i[0])])
    elif i[1] == 3:
        loop_new[1].append(loop[int(i[0])])
    elif i[1] == 7:
        loop_new[2].append(loop[int(i[0])])
    elif i[1] == 6:
        loop_new[3].append(loop[int(i[0])])
    elif i[1] == 0:
        loop_new[4].append(loop[int(i[0])])
    elif i[1] == 5:
        loop_new[5].append(loop[int(i[0])])
    elif i[1] == 2:
        loop_new[6].append(loop[int(i[0])])
    elif i[1] == 4:
        loop_new[7].append(loop[int(i[0])])

for k , v in loop_new.items():
    v = np.array(v , dtype = loop_type)
    loop_new[k] = v


#---------------------------------------------gene_expression_log2(FPKM)----------------------------------------------------------
fpkm_type = ({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_1' , 'CCS_2' , 'CCS_3' , 'fESC_1' , 'fESC_2' , 'fESC_3' , 'NT5_1'	, 'NT5_2' ,	'NT5_3' , 'NT5_4' , 'NT6_1' , 'NT6_2' , 'NT6_3'],
             'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float ]})

gene_type = ({'names':['CCS' , 'NT5' , 'NT6' , 'fESC'],
             'formats':[np.float , np.float , np.float , np.float]})    
gene_fpkm = np.loadtxt('F:\\xxli\\data_analysis\\BDF1\\RNA_seq\\gene_expression\\all_genes.txt' , skiprows = 1 , dtype = fpkm_type)

gtf_1 = []
for i in gene_fpkm:
    if i['strand'] == '+':
        gtf_1.append((i['gene_id'].split(".")[0] , i[1] , i[2] , i[3] , i[4] , i[5] , i[6] , i[7] , i[8] , i[9] , i[10] , i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17] , i[18]))
    elif i['strand'] == '-':
        gtf_1.append((i['gene_id'].split(".")[0] , i[1] , i[2] , i[3] , i[5] , i[4] , i[6] , i[7] , i[8] , i[9] , i[10] , i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17] , i[18]))
    else:
        print 1
        continue


gene_fpkm = np.array(gtf_1 , dtype = fpkm_type)

classify_gene = {}

for c in loop_new:
    classify_gene[c] = []
    for i in loop_new[c]:
        chro = 'chr' + i['chr']
        start = i['start']
        end = i['end']
        lib = gene_fpkm[gene_fpkm['chr'] == chro]
        mask = (start <= lib['end']) & (end >= lib['start'])
        overlap = lib[mask]
        if overlap.size != 0:
            for j in overlap:
                dis = min(abs(start - j['start']) , abs(end - j['end']))
                ccs = np.log2((j['CCS_1'] + j['CCS_2'] + j['CCS_3']) / 3 + 0.0001)
                nt5 = np.log2((j['NT5_1'] + j['NT5_2'] + j['NT5_3'] + j['NT5_4']) / 4 + 0.0001)
                nt6 = np.log2((j['NT6_1'] + j['NT6_2'] + j['NT6_3']) / 3 + 0.0001)
                fesc = np.log2((j['fESC_1'] + j['fESC_2'] + j['fESC_3']) / 3 + 0.0001)
                if ((ccs > 0) or (nt5 > 0) or (nt6 > 0) or (fesc > 0)) and dis <= 100000:
                    classify_gene[c].append((ccs , nt5 , nt6 , fesc))
        else:
            continue


for k , v in classify_gene.items():
    v = np.array(v , dtype = gene_type)
    classify_gene[k] = v


for c in range(8):
    left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot([classify_gene[c]['CCS'] , classify_gene[c]['NT5'] , classify_gene[c]['NT6'] , classify_gene[c]['fESC']])
    x = ['CCS','NT5','NT6','fESC']

#    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 25)
    ax.set_ylabel('log2(FPKM)' , fontsize = 25)
    ax.set_ylim((-10 , 10))
    plt.title('Gene expression(log2(FPKM)) Classify_' + str(c),fontsize = 25)
    plt.savefig('D:\\workspace\\loop\\scatter\\loop_strength_new\\loop_strength_20K_union\\gene_classify_100K\\gene_classify_' + str(c) + '.png')
    
#------------------------------------------------ATAC_Peaks------------------------------------------------------------
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
peak_type = ({'names':['chr' , 'start' , 'end' , 'score'],
             'formats':['S8' , np.int , np.int , np.float]})
gene_type = ({'names':['CCS' , 'NT5' , 'NT6' , 'fESC'],
             'formats':[np.float , np.float , np.float , np.float]})    

peaks = {}    
for c in cell:
    peak = np.loadtxt('D:\\workspace\\IDR\\selected_peaks\\idr_select_peaks_chr\\' + c + '_IDR_Selected.narrowPeak' ,
                      skiprows = 0 , usecols = (0 , 1 , 2 , 4 ) , dtype = peak_type )
    peaks[c] = peak


union = union_peaks(peaks)

classify_peaks = {}

for cla in loop_new:
    classify_peaks[cla] = []
    for i in loop_new[cla]:
        chro = i['chr']
        start = i['start']
        end = i['end']
        lib = union[union['chr'] == chro]
        mask = (start <= lib['end']) & (end >= lib['start'])
        overlap = lib[mask]
        if overlap.size != 0:
            dis_score = {}
            for j in overlap:
                dis = min(abs(start - j['start']) , abs(end - j['end']))
                if ((j['CCS'] > 1) or (j['NT5'] > 1) or (j['NT6'] > 1) or (j['fESC'] > 1)):
                    classify_peaks[cla].append((np.log2(j['CCS'] + 0.001) , np.log2(j['NT5'] + 0.001) , np.log2(j['NT6'] + 0.001) , np.log2(j['fESC'] + 0.001)))

#                    dis_score[dis] = (np.log2(j['CCS'] + 0.001) , np.log2(j['NT5'] + 0.001) , np.log2(j['NT6'] + 0.001) , np.log2(j['fESC'] + 0.001))
#            try:
#                classify_peaks[cla].append(dis_score[min(dis_score.keys())])
#            except:
#                pass
#        else:
#            pass

for k , v in classify_peaks.items():
    v = np.array(v , dtype = gene_type)
    classify_peaks[k] = v


for c in range(8):
    left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot([classify_peaks[c]['CCS'] , classify_peaks[c]['NT5'] , classify_peaks[c]['NT6'] , classify_peaks[c]['fESC']])
    x = ['CCS','NT5','NT6','fESC']

#    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 25)
    ax.set_ylabel('log2(FPKM)' , fontsize = 25)
    ax.set_ylim((-12 , 15))
    plt.title('ATAC Peak(log2(Peak_score)) Classify_' + str(c),fontsize = 25)
    plt.savefig('D:\\workspace\\loop\\scatter\\loop_strength_new\\loop_strength_20K_union\\atac_classify_all\\ATAC_peaks_classify_' + str(c) + '.png')