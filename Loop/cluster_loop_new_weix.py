# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 21:22:50 2021

@author: xxli
"""

from __future__ import  division
import numpy as np
import  math
import time
times = time.time()


def Correct_VC(X, alpha):
    x = np.array(X,float)
    s1 = np.sum(x, axis = 1)
    s1 = s1 ** alpha
#    s /= np.mean(s[s!=0])
    s1[s1 == 0] = 1
    s2 = np.sum(x, axis = 0)
    s2 = s2 ** alpha
#    s2 /= np.mean(s2[s2 != 0])
    s2[s2 == 0] = 1
    return x / (s2[None, :] * s1[:, None])


def Normal_VC_Correct(NPZ):
    """
    """
    Raw_Lib = np.load(NPZ)
    Nor_Lib = {}
    for c in Raw_Lib.keys():
        Raw_M = Raw_Lib[c]
        Nor_M = Correct_VC(Raw_M, 2/3)
        
        #Recall
        Re_factor = Raw_M.mean() / Nor_M.mean()
        Nor_M = Nor_M * Re_factor
        Nor_Lib[c] = Nor_M
    
    return Nor_Lib




def center_sites(lists):
    sum_x = 0
    sum_y = 0
    for x in lists:
        sum_x += x[1]
        sum_y += x[2]
    n = len(lists)
    return [float(sum_x)/n, float(sum_y)/n]

def distance(sites_1,sites_2):
    dis = math.sqrt((sites_1[0]-sites_2[1])**2 + (sites_1[1]-sites_2[2])**2)
    return dis

def clustering(loop_sites, dis):
    classes = []
    for i in chrom:
        c_loops = sorted(loop_sites[loop_sites['chr'] == i], key = lambda x:x[1])
        while True :
            c_classs = []
            c_classs.append(c_loops[0])
            c_loops.remove(c_loops[0])
            center = center_sites(c_classs)
            for loop in c_loops:     
                 if distance(center, loop) <= dis:
                      c_classs.append(loop)
                      center = center_sites(c_classs)
                      c_loops.remove(loop)
            classes.append(c_classs)
            if len(c_loops) == 0:
                break
    return classes





def filtering(clss):
    loop = []
    for outs in clss:
        lens = len(outs)
        if  len(outs) >= 2 :
            tmp = sorted(outs, key=lambda x:x[-1])
            n = [] ; IF = []
            for i in range(len(tmp)):
                if tmp[i][-1] == tmp[0][-1]:
                    n.append(i) 
                    IF.append(tmp[i][-2])
            index = np.argmax(IF)
            x = n[index]
            # print (outs)
            tmps = (tmp[x][0] , tmp[x][1] , tmp[x][2] , tmp[x][3] , tmp[x][4] / (10**lens))
            loop.append(tmps)
        else :
            loop.append(outs[0])
    return np.array(loop, dtype = loopType)    



    
# -----------------------Greedy for clustering loops-------------------------------
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']
#load loop Data
cell = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'NTs' , 'fESC']
res = '20K'
R = 20000
initial_distance = 2 * R * math.sqrt(2) + 1000
loopType = np.dtype({'names':['chr' , 'S1' , 'E1' , 'IF' , 'Q-value'],
                     'formats':['S8', np.int, np.int , np.float , np.float]})

#---------------------------------Files--------------------------
CCS_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/CCS_matrix/CCS_20K.npz')
NT5_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/NT5_matrix/NT5_20K.npz')
NT6_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/NT6_matrix/NT6_20K.npz')
F35_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/F35_matrix/F35_20K.npz')
F40_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/F40_matrix/F40_20K.npz')
NTs_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/NTs_matrix/NTs_20K.npz')
fESC_Lib = Normal_VC_Correct('/public/home/lixinxin/data/BDF1/HiC/HapHiC_workspace/Cooler_to_matrix/raw/fESC_matrix/fESC_20K.npz')





HiC_Data = {'CCS':CCS_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib,
            'F35':F35_Lib,
            'F40':F40_Lib,
            'NTs':NTs_Lib,
            'fESC':fESC_Lib}







for c in cell:
    Lib = HiC_Data[c]
    loopSource = '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/' + c + '_loops_' + res + '.txt'
    loop_out = '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_' + c + '_loops_'  + res + '_3_0.0001_new.txt'

    Loops = np.loadtxt(loopSource, dtype = loopType, usecols = [0 , 1 , 2 , 3 , -1], skiprows = 1)
    mask = Loops['Q-value'] <= 0.01
    Loops = Loops[mask]
    print len(Loops)
    
    loops = []
    for i in Loops:
        chro = i['chr']
        matrix = Lib[chro]
        startHiC = i['S1'] // R
        endHiC = i['E1'] // R
        loops.append((i['chr'] , i['S1'] , i['E1'] , matrix[startHiC , endHiC] , i['Q-value']))
    loops = np.array(loops , dtype = loopType)

        
            
    cluster_1 = clustering(loops, initial_distance)
    #loops2 = []
    #for outs in cluster_1:
    #    lens = len(outs)
    #    if lens >= 2 :
    #        c_center = center_sites(outs)
    #        cuts = initial_distance
    #        for i in range(len(outs)):
    #            dis = distance(c_center, outs[i])
    #            if dis < cuts:
    #                cuts = dis
    #                n = i
    #                outs[n][3] = outs[n][3]/(lens**3)
    #        loops2.append(outs[n])   
    #    elif len(outs) == 1 and outs[0]['Q-value'] <= 0.005:
    #        loops2.append(outs[0])
    #        
    #loops2 = np.array(loops2, dtype = loopType)
    loops2 = filtering(cluster_1)
    cluster_2 = clustering(loops2, initial_distance)
    
    loops3 = filtering(cluster_2)
    cluster_3 = clustering(loops3, initial_distance)  
    
    loops_out = filtering(cluster_3)               
                    
    file_out = open(loop_out,'w')
    file_out.writelines('chr\tstart\tend\tQ-value\n')
    for outs in loops_out:
         if ((outs[4] < 0.0001) and (outs[3] > 50)) or (outs[4] == 0):
             strs = '\t'.join(map(str,outs)) + '\n'
             file_out.write(strs)
    file_out.close()    
    print len(loops)
    
    print (time.time() - times  )
    
    
    
    
    
    
for c in ['NTs' , 'fESC']:
    Lib = HiC_Data[c]
    loopSource = '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/' + c + '_loops_' + res + '.txt'
    loop_out = '/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/Clustered_' + c + '_loops_'  + res + '_3_0.000001.txt'

    Loops = np.loadtxt(loopSource, dtype = loopType, usecols = [0 , 1 , 2 , 3 , -1], skiprows = 1)
    mask = Loops['Q-value'] <= 0.01
    Loops = Loops[mask]
    print len(Loops)
    
    loops = []
    for i in Loops:
        chro = i['chr']
        matrix = Lib[chro]
        startHiC = i['S1'] // R
        endHiC = i['E1'] // R
        loops.append((i['chr'] , i['S1'] , i['E1'] , matrix[startHiC , endHiC] , i['Q-value']))
    loops = np.array(loops , dtype = loopType)

        
            
    cluster_1 = clustering(loops, initial_distance)
    #loops2 = []
    #for outs in cluster_1:
    #    lens = len(outs)
    #    if lens >= 2 :
    #        c_center = center_sites(outs)
    #        cuts = initial_distance
    #        for i in range(len(outs)):
    #            dis = distance(c_center, outs[i])
    #            if dis < cuts:
    #                cuts = dis
    #                n = i
    #                outs[n][3] = outs[n][3]/(lens**3)
    #        loops2.append(outs[n])   
    #    elif len(outs) == 1 and outs[0]['Q-value'] <= 0.005:
    #        loops2.append(outs[0])
    #        
    #loops2 = np.array(loops2, dtype = loopType)
    loops2 = filtering(cluster_1)
    cluster_2 = clustering(loops2, initial_distance)
    
    loops3 = filtering(cluster_2)
    cluster_3 = clustering(loops3, initial_distance)  
    
    loops_out = filtering(cluster_3)               
                    
    file_out = open(loop_out,'w')
    file_out.writelines('chr\tstart\tend\tQ-value\n')
    for outs in loops_out:
         if ((outs[4] < 0.000001) and (outs[3] > 60)) or (outs[4] == 0):
             strs = '\t'.join(map(str,outs)) + '\n'
             file_out.write(strs)
    file_out.close()    
    print len(loops)
    
    print (time.time() - times  )
    
    
    
    