# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 22:36:58 2020

@author: han-luo
"""

import matplotlib
import matplotlib.pyplot as plt




def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2.- 0.2, 1.03*height, '%s' % int(height))


name_list = ['CCS', 'NT5', 'NT6', 'fESC']
num_list = [2724, 802, 946, 604]


autolabel(plt.bar(range(len(num_list)), num_list, color='lightseagreen', tick_label=name_list))
plt.show()