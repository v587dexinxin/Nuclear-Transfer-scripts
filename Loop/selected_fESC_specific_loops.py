# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 03:18:22 2021

@author: xxli
"""

    
def apa_submatrix(M, pos, w=5):
    
    Len = M.shape[0]

    apa = []
    for i, j in pos:
        if (i-w>=0) and (i+w+1<=Len) and (j-w>=0) and (j+w+1<=Len):
            tmp = M[i-w:i+w+1, j-w:j+w+1]
            if tmp.mean()==0:
                continue
            mask = np.isnan(tmp)
            if mask.sum() > 0:
                continue
            tmp = tmp / tmp.mean()
            apa.append(tmp)
    
    return apa

def apa_analysis(apa, w=5, cw=3):
    
    avg = apa.mean(axis=0)
    lowerpart = avg[-cw:,:cw]
    upperpart = avg[:cw,-cw:]
    maxi = upperpart.mean() * 5
    ## APA score
    score = avg[w,w] / lowerpart.mean()
    ## z-score
    z = (avg[w,w] - lowerpart.mean()) / lowerpart.std()
    p = 1 - ndtr(z)
    
    return avg, score, z, p, maxi

data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S64' , np.int , np.int]})


chros = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']



NT_fESC_speci_Loops = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/fESC_Speci_Loops_NT_fESC_merged_20K.txt' , 
                             dtype = data_type , usecols = (0 , 1 , 2))




out = open('/public/home/lixinxin/data/BDF1/HiC/Loop/Loop_new/HiCCUPs/Clustered_new/classify_loops_weix/Selected_NT_fESC_specific_loops_20K.txt' , 'w')

n = 0
                  
for g in chros:
    peaks = NT_fESC_speci_Loops[NT_fESC_speci_Loops['chr'] == g]
    
    M = HiC_Data['CCS'][g]
    for p in peaks:
        pos = []
        apa = []  
        x, y = p[1], p[2]
        if abs(y-x) < 15 * res:
            continue
        s_l = range(p[1]//res, int(np.ceil((p[1]+20000)/float(res))))
        e_l = range(p[2]//res, int(np.ceil((p[2]+20000)/float(res))))
        si, ei = None, None
        for st in s_l:
            for et in e_l:
                if (st < M.shape[0]) and (et < M.shape[0]):
                    if si is None:
                        si, ei = st, et
                    else:
                        if M[st,et]>M[si,ei]:
                            si,ei = st,et
        if not si is None:
            if si < ei:
                pos.append((si,ei))
            else:
                pos.append((ei,si))
        tmp = apa_submatrix(M, pos)
        if len(tmp) == 0:
            print p
            continue
        apa.extend(tmp)
        apa = np.r_[apa]
        avg,score,z,p_value,maxi = apa_analysis(apa)
        if score > 5:
            n += 1
        else:
               out.writelines('\t'.join([p[0] , str(p[1]) , str(p[2])]) + '\n')
out.close()

