#!/usr/bin/env python

##this script we will generate black list from blasting results
import re
import glob
import sys
import subprocess
import os

input_blast_opt_fl = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/00_make_black_list_042821/07_make_black_blast_mtpt_genome_050321/01_blast_mtpt_genome_050321/output_dir/opt_blast.txt

input_output_dir = sys.argv[2]

input_organelle_str = sys.argv[3]

thr = 0.8

def make_black (input_blast_opt_fl,input_output_dir,input_organelle_str):

    store_final_black_line_list = []
    total_len = 0
    with open (input_blast_opt_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            subject_chr = col[1]
            iden = col[2]
            subject_loc = subject_chr + '\t' + col[8] + '\t' + col[9]
            len = int(col[4])

            if subject_chr not in input_organelle_str.split(','):

            #if subject_chr != 'ChrPt' and subject_chr != 'ChrMt':

                if float(iden) > thr:
                    store_final_black_line_list.append(subject_loc)

                    total_len = total_len + len

    with open (input_output_dir + '/opt_mtpt_black.bed','w+') as opt:
        for eachline in store_final_black_line_list:
            opt.write(eachline + '\n')

    #with open (input_output_dir + '/opt_total_black_len.txt','w+') as opt:
    #    opt.write(str(total_len))

make_black (input_blast_opt_fl,input_output_dir,input_organelle_str)



