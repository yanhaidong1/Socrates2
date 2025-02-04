#!/usr/bin/env python

##this script is to run the fimo using the fimo based on the ACRs

import sys
import re
import glob
import os
import subprocess
from multiprocessing import Pool
#from Bio import SeqIO

#input_genome_fasta_fl = sys.argv[1] ##contains the fasta information for each chr
#input_motif_fl = sys.argv[2]
#input_acr_fl = sys.argv[3]


#chr_list_fl = sys.argv[3] ##this file contains the chrom information that could be divided into multiple chromosome
#output_dir = sys.argv[4] ##create genome

#extend_size = sys.argv[5]
#process_num = sys.argv[6]
#pval_cutoff = sys.argv[7]


def multi_run(args):
    return analyze_acr_motif (*args)

def analyze_acr_motif(ipt_acr_fl,ipt_genome_fasta_fl,ipt_motif_fl,opt_dir,extend_size):

    store_final_line_list = []
    with open (ipt_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            extend_st = int(col[1]) - int(extend_size)
            extend_end = int(col[2]) + int(extend_size)

            if extend_st < 0:
                final_extend_st = 0
            else:
                final_extend_st = extend_st

            final_line = col[0] + '\t' + str(final_extend_st) + '\t' + str(extend_end) + '\t' + col[0] + '_' + col[1] + '_' + col[2]
            store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_final_acr_extended.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_final_acr_extended.txt > ' + \
          opt_dir + '/temp_final_acr_extended_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##extract fasta of the final acr
    cmd = 'bedtools getfasta -fi ' + ipt_genome_fasta_fl + \
          ' -bed ' +  opt_dir + '/temp_final_acr_extended_sorted.txt' + \
          ' -name > ' +  opt_dir + '/temp_final_acr_extended.fasta'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##updating 010324
    ##we will update the fasta file to remove the location of the title
    store_final_line_list = []
    with open (opt_dir + '/temp_final_acr_extended.fasta','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if eachline.startswith('>'):

                if re.match('(.+)::.+',eachline):
                    mt = re.match('(.+)::.+',eachline)
                    final_line = mt.group(1)
                else:
                    final_line = eachline

                store_final_line_list.append(final_line)
            else:
                store_final_line_list.append(eachline)

    with open (opt_dir + '/temp_final_acr_extended_modi.fasta','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##run the fimo
    cmd = 'fimo -max-stored-scores 10000000 --oc ' + opt_dir + ' ' + ipt_motif_fl + ' ' + opt_dir + '/temp_final_acr_extended_modi.fasta'
    print(cmd)
    subprocess.call(cmd, shell=True)


def pipeline_analyze_acr_motif (input_acr_fl,input_genome_fasta_fl,input_motif_fl,output_dir,
                                process_num,extend_size):

    temp_all_output_dir = output_dir + '/temp_all_output_dir'
    if not os.path.exists(temp_all_output_dir):
        os.makedirs(temp_all_output_dir)

    temp_split_acr_dir = output_dir + '/temp_split_acr_dir'
    if not os.path.exists(temp_split_acr_dir):
        os.makedirs(temp_split_acr_dir)

    for x in range(0, int(process_num)):
        dir_code = x + 1
        temp_output_dir = temp_all_output_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)


    ##split the acr file
    cmd = 'split -n l/' + str(process_num) + ' ' + input_acr_fl + ' ' + temp_split_acr_dir + '/temp_split_'
    print(cmd)
    subprocess.call(cmd,shell=True)

    temp_split_file_list = glob.glob(temp_split_acr_dir + '/temp_split*')
    temp_all_output_dir_list = glob.glob(temp_all_output_dir + '/*')

    pool = Pool(int(process_num))
    run_list = []
    for x in range(0, int(process_num)):
        each_func_argument = (temp_split_file_list[x],
                              input_genome_fasta_fl,
                              input_motif_fl,
                              temp_all_output_dir_list[x],
                              extend_size)
        run_list.append(each_func_argument)
    pool.map(multi_run, run_list)


##later we need to generate a function to combine all the information together
def wrap_all_results (output_dir,pval_cutoff):

    ##since the gff file is so large we need to transfer to bed first and combine all together
    ##generate a dir to store the bed
    temp_store_gff2bed_dir = output_dir + '/temp_store_gff2bed_dir'
    if not os.path.exists(temp_store_gff2bed_dir):
        os.makedirs(temp_store_gff2bed_dir)

    store_all_bed_fl_list = []
    target_output_dir = output_dir + '/temp_all_output_dir'
    all_opt_dir_list = glob.glob(target_output_dir + '/*')
    for eachdir in all_opt_dir_list:

        mt = re.match('.+/(.+)',eachdir)
        dirnm = mt.group(1)
        mt = re.match('temp_output_(\d+)dir',dirnm)
        rep_nm = mt.group(1)

        gff_fl = eachdir + '/fimo.gff'

        gff_fl_o = open(gff_fl,'r')
        bed_fl_w = open(temp_store_gff2bed_dir + '/temp_bed_' + rep_nm + '.bed','w')

        for eachline in gff_fl_o:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if not eachline.startswith('#'):
                acr_loc = col[0]
                st = col[3]
                ed = col[4]

                annot_col = col[8].split(';')
                #qval = '1'
                #motifID = ''
                #seq = ''

                store_annot_dic = {}
                for eachcol in annot_col:
                    if '=' in eachcol:
                        annot_col2 = eachcol.split('=')
                        store_annot_dic[annot_col2[0]] = annot_col2[1]

                if 'Name' in store_annot_dic:
                    motifID = store_annot_dic['Name']
                    mt = re.match('(.+)_.+_.+_.+',motifID)
                    final_motifID = mt.group(1)

                else:
                    final_motifID = 'unknown'
                if 'pvalue' in store_annot_dic:
                    pval = store_annot_dic['pvalue']
                else:
                    pval = 'unknown'
                if 'sequence' in store_annot_dic:
                    seq = store_annot_dic['sequence']
                else:
                    seq = 'unknown'

                final_line = acr_loc + '\t' + st + '\t' + ed + '\t' + final_motifID + '\t' + pval + '\t' + col[6] + '\t' + seq
                bed_fl_w.write(final_line + '\n')

        gff_fl_o.close()
        bed_fl_w.close()

        store_all_bed_fl_list.append(temp_store_gff2bed_dir + '/temp_bed_' + rep_nm + '.bed')

    ##combine all the bed together
    bed_fl_str = ' '.join(store_all_bed_fl_list)
    cmd = 'cat ' + bed_fl_str + ' > ' + output_dir + '/opt_final_motif.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##sort the bed file
    cmd = 'sort -k1,1V -k2,2n ' + output_dir + '/opt_final_motif.bed' + ' > ' + output_dir + '/opt_final_motif_sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)


    motif_fl = open (output_dir + '/opt_final_motif_sorted.bed','r')
    rm_dup_motif_fl = open(output_dir + '/opt_final_motif_rmdup_flt_sorted.bed','w')

    store_line_dic = {}
    for eachline in motif_fl:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        if float(col[4]) < float(pval_cutoff):
            motif_loc = col[0] + '_' + col[1] + '_' + col[2]

            if motif_loc not in store_line_dic:
                rm_dup_motif_fl.write(eachline + '\n')

            store_line_dic[motif_loc] = 1

    motif_fl.close()
    rm_dup_motif_fl.close()

    ##build a ACR motif file
    store_final_line_list = []
    with open (output_dir + '/opt_final_motif_rmdup_flt_sorted.bed','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr_location = '\t'.join(col[0].split('_'))
            chrnm = col[0].split('_')[0]

            final_line = acr_location + '\t' + 'na' + '\t' + 'na' + '\t' + chrnm + '\t' + col[1] + '\t' + col[2] + '\t' + \
                         col[3] + '\t' + col[4] + '\t' + col[5] + '\t' + col[6]
            store_final_line_list.append(final_line)

    with open (output_dir + '/opt_final_ACR_motif.bed', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



#pipeline_analyze_acr_motif (input_acr_fl,input_genome_fasta_fl,input_motif_fl,output_dir,
#                                process_num,extend_size)
#wrap_all_results (output_dir)








