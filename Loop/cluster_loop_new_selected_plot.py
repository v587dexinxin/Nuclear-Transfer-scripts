# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:58:21 2021

@author: xxli
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 21:17:44 2019

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
        lens = len(outs)
        if  len(outs) >= 2 :
            tmp = sorted(outs, key=lambda x:x[4])
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

for c in cell:
    loopSource = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS\\' + c + '_loops_' + res + '.txt'
    loop_out = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS_new\\Clustered_loops_selected_plot\\Clustered_' + c + '_loops_'  + res + '_3_0.0001.txt'
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
    cluster_2 = clustering(loops2, initial_distance)
    
    loops3 = filtering(cluster_2)
    cluster_3 = clustering(loops3, initial_distance)  
    
    loops_out = filtering(cluster_3)               
                    
    file_out = open(loop_out,'w')
    file_out.writelines('chr\tstart\tend\tQ-value\n')
    for outs in loops_out:
         if ((outs[4] < 0.0001) and (outs[3] > 47)):
             strs = '\t'.join(map(str,outs)) + '\n'
             file_out.write(strs)
    file_out.close()    
    
    print (time.time() - times  )