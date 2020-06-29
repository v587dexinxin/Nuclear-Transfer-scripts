# -*- coding: utf-8 -*-
"""
Created on Sat Mar 03 11:45:49 2018

@author: xxli
"""

f1 = open('C:\\Users\\xxli\\Desktop\\mm10.txt','r')
f2 = open('C:\\Users\\xxli\\Desktop\\mm10_10K.bed','w')
result = {}
R = 10000
for i in f1:
    i = i.split()
    result[i[0]] = []
    binnum = int(i[1])/R
    for j in range(binnum):
        result[i[0]].append([j*R,(j+1)*R])
        f2.writelines(i[0] + '\t' + str(j*R) + '\t' + str((j+1)*R) + '\n')
f1.close()
f2.close()