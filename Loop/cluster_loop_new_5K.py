# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:51:21 2021

@author: xxli
"""

from __future__ import  division
import numpy as np
import  math
import time
times = time.time()


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
        lens = len(outs) - 1
        if  len(outs) >= 2 :
            tmp = sorted(outs, key=lambda x:x[4])
            # print (outs)
            tmp = (tmp[0][0] , tmp[0][1] , tmp[0][2] , tmp[0][3] , tmp[0][4] / (10**lens))
            loop.append(tmp)
        else :
            loop.append(outs[0])
    return np.array(loop, dtype = loopType)    



    
# -----------------------Greedy for clustering loops-------------------------------
chrom = ['4']
#load loop Data
cell = ['CCS' , 'NTs' , 'fESC']
res = '5K'
R = 5000

for c in cell:
    loopSource = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS_new\\' + c + '_loops_' + res + '.txt'
    loop_out = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS_new\\Clustered_' + c + '_loops_'  + res + '_3_0.01.txt'
    loopType = np.dtype({'names':['chr' , 'S1' , 'E1' , 'IF' , 'Q-value'],
                         'formats':['U8', np.int, np.int , np.float , np.float]})
    loops = np.loadtxt(loopSource, dtype = loopType, usecols = [0 , 1 , 2 , 3 , -1], skiprows = 1)
    mask = loops['Q-value'] <= 0.01
    loop_hq = loops[mask].copy()
    initial_distance = 2 * R * math.sqrt(2) + 1000
    

        
            
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
    cluster_2 = clustering(loops2, initial_distance * 2)
    
    loops3 = filtering(cluster_2)
    cluster_3 = clustering(loops3, initial_distance * 2)  
    
    loops_out = filtering(cluster_3)               
                    
    file_out = open(loop_out,'w')
    file_out.writelines('chr\tstart\tend\tQ-value\n')
    for outs in loops_out:
         if (outs[4] < 0.01) and (outs[3] >= 50):
             strs = '\t'.join(map(str,outs)) + '\n'
             file_out.write(strs)
    file_out.close()    
    
    print (time.time() - times)
    
    
    
    
    
    
    
    
    
    
    
    
    