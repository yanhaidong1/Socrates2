#!/usr/bin/env python

##this script we will generate a black list from the repeatmasker results

import re
import glob
import sys
import subprocess
import os


input_repeat_masker_fl = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/Osativa323v7.fa.out.gff
#/scratch/hy17471/rice_altlas_scATAC_seq_042021/09_3_crossSpecies_analysis_trajectory_cellAlign_031723/00_analysis_on_Arabidopsis_041823/raw/repeatmasker_run_042023/Athaliana_167_TAIR10.fa.out.gff

input_output_dir = sys.argv[2]

def make_black_list (input_repeat_masker_fl,input_output_dir):

    store_final_line_list = []
    length = 0
    with open (input_repeat_masker_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()
                tenm = col[9]

                if re.match('\"Motif:\((.+)\)n"',tenm):
                    length =  length + (int(col[4]) - int(col[3]))
                    final_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(final_line)

    #with open (input_output_dir + '/opt_repeat_length.txt','w+') as opt:
    #    opt.write(str(length))

    with open (input_output_dir + '/opt_repeat_black.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

make_black_list (input_repeat_masker_fl,input_output_dir)
