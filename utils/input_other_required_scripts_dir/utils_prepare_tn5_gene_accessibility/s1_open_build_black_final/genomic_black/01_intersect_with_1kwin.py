#!/usr/bin/env python

##this script we will generate 1kb window of genome and intersect it with the tn5
import re
import glob
import sys
import subprocess
import os


input_tn5_bed_fl = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/00_make_black_list_042821/04_change_bam_to_bed_042921/output_dir_atac/atac_input/opt_tn5_sorted_unique.bed

input_genome_fai_fl = sys.argv[2]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/Osativa323v7.fa.fai

input_nb_fastSparsetn5_pl = sys.argv[3]
##fastSparse.nonbinary.peak.pl

input_output_dir = sys.argv[4]


def make_win_intersect (input_tn5_bed_fl,input_genome_fai_fl,input_nb_fastSparsetn5_pl,
                        input_output_dir):

    ##generate 1kbwindows for the genome.fa.fai
    cmd = 'bedtools makewindows -g ' + input_genome_fai_fl + ' -w 1000 > ' + input_output_dir + '/1kbwindows.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -a ' + input_tn5_bed_fl + \
          ' -b ' + input_output_dir + '/1kbwindows.bed' + ' -wa -wb' + \
          ' -g ' + input_genome_fai_fl + \
          ' -sorted > ' + input_output_dir + '/temp_intersect.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'perl ' + input_nb_fastSparsetn5_pl + \
          ' ' + input_output_dir + '/temp_intersect.txt' + ' > ' + input_output_dir + '/opt_nb_win_read.sparse'
    print(cmd)
    subprocess.call(cmd, shell=True)

make_win_intersect (input_tn5_bed_fl,input_genome_fai_fl,input_nb_fastSparsetn5_pl,
                        input_output_dir)