#!/usr/bin/env python

##this script we will identify the black regions from the tn5 number in the atac and chip

import re
import glob
import sys
import subprocess
import os


input_atac_fl = sys.argv[1]

#input_chip_fl = sys.argv[2]

input_output_dir = sys.argv[2]

def identify_black_region (input_atac_orchip_fl):

    ##calculate the average of tn5 for all window
    store_win_read_num_dic = {}
    with open (input_atac_orchip_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            win = col[0]
            read_num = col[2]

            if win in store_win_read_num_dic:
                store_win_read_num_dic[win] += int(read_num)
            else:
                store_win_read_num_dic[win] = int(read_num)

    avg_num_win = 0
    win = 0
    for eachwin in store_win_read_num_dic:
        win_num = store_win_read_num_dic[eachwin]
        avg_num_win = avg_num_win + win_num
        win += 1

    print(avg_num_win)
    print(win)

    f4_avg = 4*(avg_num_win/win)


    store_win_over_avg_line_list = []
    store_all_line_list = []
    for eachwin in store_win_read_num_dic:
        if store_win_read_num_dic[eachwin] >= f4_avg:
            mt = re.match('(.+)_(.+)_(.+)',eachwin)
            chr_nm = mt.group(1)
            st = mt.group(2)
            ed = mt.group(3)

            final_line = chr_nm + '\t' + st + '\t' + ed
            store_win_over_avg_line_list.append(final_line)

        final_line = eachwin + '\t' + str(store_win_read_num_dic[eachwin])
        store_all_line_list.append(final_line)

    return (store_win_over_avg_line_list,store_all_line_list)

atac_line_list,store_all_line_list = identify_black_region (input_atac_fl)

with open (input_output_dir + '/opt_genomic_black.bed','w+') as opt:
    for eachline in atac_line_list:
        opt.write(eachline + '\n')

#with open (input_output_dir + '/opt_atac_all_win_num.txt','w+') as opt:
#    for eachline in store_all_line_list:
#        opt.write(eachline + '\n')


#chip_line_list,store_all_line_list = identify_black_region (input_chip_fl)

#with open (input_output_dir + '/opt_chip_black.bed','w+') as opt:
#    for eachline in chip_line_list:
#        opt.write(eachline + '\n')

#with open (input_output_dir + '/opt_chip_all_win_num.txt','w+') as opt:
#    for eachline in store_all_line_list:
#        opt.write(eachline + '\n')
