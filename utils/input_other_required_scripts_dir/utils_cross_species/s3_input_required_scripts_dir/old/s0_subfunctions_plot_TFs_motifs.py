#!/usr/bin/env python


##Here we will check the TF accessiblity

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats
from Bio import SeqIO
from statistics import mean
from operator import itemgetter



##plot motifs or TFs
def subfunction_plot_TFs (input_required_scripts_dir,ipt_TF_motif_corr_fl,ipt_gene_impute_acc_rds_fl,GAobj_rds_fl
                          ,opt_dir,s0_s5_sub_lim,s0_s5_sub_plot_topnum_pairs):

    ##opt_corresponding_motif_TF.txt is the ipt_TF_motif_corr_fl
    ##s0_s5_sub_plot_topnum_pairs is the number of
    store_TF_motif_pair_coef_dic = {}
    with open (ipt_TF_motif_corr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motif_TF = col[0] + '__' + col[1]
            coef = float(col[2])
            store_TF_motif_pair_coef_dic[motif_TF] = coef

    ##select the top for building the candidate markers
    res = dict(sorted(store_TF_motif_pair_coef_dic.items(), key=itemgetter(1), reverse=True)[:int(s0_s5_sub_plot_topnum_pairs)])

    store_target_gene_dic = {}
    for eachpair in res:
        mt = re.match('(.+)__(.+)',eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_gene_dic[gene] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' +  'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachgene in store_target_gene_dic:
        final_line = eachgene + '\t' + eachgene + '\t' + 'TF' + '\t' + 'all' + '\t' + eachgene
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_TF_markers.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    plot_TFs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_genes_112723.R'

    cmd = 'Rscript ' + plot_TFs_script + \
          ' ' + ipt_gene_impute_acc_rds_fl + \
          ' ' + GAobj_rds_fl + \
          ' ' + opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_TF_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd,shell=True)


##plot motifs or TFs
def subfunction_plot_motifs (input_required_scripts_dir,ipt_TF_motif_corr_fl,ipt_motif_impute_acc_threesparse_fl, ipt_meta_fl,
                          opt_dir,s0_s5_sub_lim,s0_s5_sub_plot_topnum_pairs):

    ##opt_corresponding_motif_TF.txt is the ipt_TF_motif_corr_fl
    ##s0_s5_sub_plot_topnum_pairs is the number of
    store_TF_motif_pair_coef_dic = {}
    with open (ipt_TF_motif_corr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motif_TF = col[0] + '__' + col[1]
            coef = float(col[2])
            store_TF_motif_pair_coef_dic[motif_TF] = coef

    ##select the top for building the candidate markers
    res = dict(sorted(store_TF_motif_pair_coef_dic.items(), key=itemgetter(1), reverse=True)[:int(s0_s5_sub_plot_topnum_pairs)])

    store_target_motif_dic = {}
    for eachpair in res:
        mt = re.match('(.+)__(.+)',eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_motif_dic[motif] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' +  'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachmotif in store_target_motif_dic:
        final_line = eachmotif + '\t' + eachmotif + '\t' + 'TF' + '\t' + 'all' + '\t' + eachmotif
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_motif_markers.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    plot_motifs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_motifs_112723.R'

    cmd = 'Rscript ' + plot_motifs_script + \
          ' ' + ipt_motif_impute_acc_threesparse_fl + \
          ' ' + ipt_meta_fl + \
          ' ' + opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_motif_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd,shell=True)











