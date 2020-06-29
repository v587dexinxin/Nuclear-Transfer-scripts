# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:03:58 2019

@author: han-luo
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
apa = []
res = 40000
for c in chros:
    peaks = P_Loops[P_Loops['chr'] == c]
    pos = []
    M = P_40K[c]
    for p in peaks:
        x, y = p[1], p[2]
        if abs(y-x) < 10 * res:
            continue
        s_l = range(p[1]//res, int(np.ceil((p[1]+40000)/float(res))))
        e_l = range(p[2]//res, int(np.ceil((p[2]+40000)/float(res))))
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
    apa.extend(tmp)
apa = np.r_[apa]
avg,score,z,p,maxi = apa_analysis(apa)
plt.imshow(avg, cmap=cmap, vmax=3, interpolation='none')
plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
plt.title('APA score = {0:.3g}, p-value = {1:.3g}'.format(score, p))
plt.colorbar()