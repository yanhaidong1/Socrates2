#!/usr/bin/env python

##this script is to prepare the data for the motif enrichment test
import re
import glob
import sys
import subprocess
import os
import random



def prepare_peak_tn5_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
                             prefix,
                             input_output_dir):

    opt_sorted_peak_bed_dir = input_output_dir + '/opt_sorted_peak_bed_dir'
    if not os.path.exists(opt_sorted_peak_bed_dir):
        os.makedirs(opt_sorted_peak_bed_dir)

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'
    if not os.path.exists(opt_peak_sparse_dir):
        os.makedirs(opt_peak_sparse_dir)

    peak_bed_fl_path = input_peak_fl

    # peak_bed_fl_path = eachdir + '/' + dir_nm + '.unique500bpPeaks.bed'
    sorted_tn5_bed_fl_path = input_tn5_bed_fl

    cmd = 'sort -k1,1V -k2,2n ' + peak_bed_fl_path + ' > ' + opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'
    subprocess.call(cmd, shell=True)
    sorted_peak_bed_fl_path = opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'

    cmd = 'bedtools intersect -a ' + sorted_tn5_bed_fl_path + \
          ' -b ' + sorted_peak_bed_fl_path + ' -wa -wb' + \
          ' | perl ' + input_fastSparsetn5_pl + \
          ' - > ' + opt_peak_sparse_dir + '/opt_peak_' + prefix + '.sparse'
    print(cmd)
    subprocess.call(cmd, shell=True)

    #cmd = 'bedtools intersect -a ' + sorted_tn5_bed_fl_path + \
    #      ' -b ' + sorted_peak_bed_fl_path + ' -wa -wb' + \
    #      ' -g ' + input_genome_fai_fl + \
    #      ' -sorted | perl ' + input_fastSparsetn5_pl + \
    #      ' - > ' + opt_peak_sparse_dir + '/opt_peak_' + prefix + '.sparse'
    #print(cmd)
    #subprocess.call(cmd, shell=True)

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'

    opt_peak_sparse_sorted_dir = input_output_dir + '/opt_peak_sparse_sorted_dir'
    if not os.path.exists(opt_peak_sparse_sorted_dir):
        os.makedirs(opt_peak_sparse_sorted_dir)

    allsparse_fl_list = glob.glob(opt_peak_sparse_dir + '/*')
    for eachfl in allsparse_fl_list:
        mt = re.match('.+/(.+)\.sparse', eachfl)
        flnm = mt.group(1)

        cmd = 'sort -k1,1V -k2,2n ' + eachfl + ' > ' + opt_peak_sparse_sorted_dir + '/' + flnm + '_sorted.sparse'
        print(cmd)
        subprocess.call(cmd, shell=True)


def build_peak_cell_in_meta (input_meta_fl,input_peak_cell_sparse_fl,input_output_dir,prefix_meta):

    store_cell_peak_num_dic = {}
    with open (input_peak_cell_sparse_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            peaknm = col[0]
            cell_nm = col[1]

            if peaknm.startswith('chr'):
                mt = re.match('chr(.+)', peaknm)
                peak_nm = mt.group(1)
            else:
                peak_nm = peaknm

            if cell_nm in store_cell_peak_num_dic:

                if peak_nm in store_cell_peak_num_dic[cell_nm]:
                    store_cell_peak_num_dic[cell_nm][peak_nm] += 1
                else:
                    store_cell_peak_num_dic[cell_nm][peak_nm] = 1

            else:
                store_cell_peak_num_dic[cell_nm] = {}
                store_cell_peak_num_dic[cell_nm][peak_nm] = 1

    store_final_line_list = []
    count = 0
    with open (input_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1
            if count != 1:
                cellnm = col[0]

                if cellnm in store_cell_peak_num_dic:
                    peak_dic = store_cell_peak_num_dic[cellnm]
                    ACR_num = len(list(peak_dic.keys()))
                else:
                    ACR_num = 0

                final_line = eachline + '\t' + str(ACR_num)
                store_final_line_list.append(final_line)

            else:
                final_line = eachline + '\t' + 'ACRnum'
                store_final_line_list.append(final_line)

    with open (input_output_dir + '/opt_' + prefix_meta + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

##upating 010225
def downsampling (cutoff_cellnm,input_output_dir,prefix_meta,target_colnm):

    new_meta_fl = input_output_dir + '/opt_' + prefix_meta + '.txt'

    store_celltype_celllist_dic = {}
    count = 0
    target_index = 0
    with open (new_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1

            if count == 1:
                target_index = col.index(target_colnm)
            else:

                celltype = col[int(target_index)]

                cellnm = col[0]
                #celltype = col[24]

                if celltype in store_celltype_celllist_dic:
                    store_celltype_celllist_dic[celltype].append(cellnm)
                else:
                    store_celltype_celllist_dic[celltype] = []
                    store_celltype_celllist_dic[celltype].append(cellnm)

    store_target_cell_list = []
    for eachcelltype in store_celltype_celllist_dic:
        celllist = store_celltype_celllist_dic[eachcelltype]

        if len(celllist) > int(cutoff_cellnm):

            ##do the randomly picking
            random_picked_celllist = random.choices(celllist,k=int(cutoff_cellnm))

            for eachcell in random_picked_celllist:
                store_target_cell_list.append(eachcell)

        else:
            for eachcell in celllist:
                store_target_cell_list.append(eachcell)

    store_final_line_list = []
    count = 0
    with open (new_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                cellnm = col[0]

                if cellnm in store_target_cell_list:
                    store_final_line_list.append(eachline)

            else:
                store_final_line_list.append(eachline)

    with open (input_output_dir + '/opt_downsample_meta.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


##define a function to generate a new peak by cell sparse that will be used for generating a peak by cell matrix
##not use currenlty, as although the cell number is decreasing the aplitute does not change
def filter_peak_sparse (input_output_dir,input_peak_cell_sparse_fl):

    downsample_meta_fl = input_output_dir + '/opt_downsample_meta.txt'

    store_target_cell_dic = {}
    count = 0
    with open (downsample_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                store_target_cell_dic[col[0]] = 1


    sparse_fl = open (input_peak_cell_sparse_fl,'r')
    sparse_flt_fl = open (input_output_dir + '/opt_downsample_peak.sparse','w')

    for eachline in sparse_fl:

        eachline = eachline.strip('\n')
        col = eachline.strip().split()
        cellnm = col[1]

        peaknm = col[0]
        if peaknm.startswith('chr'):
            mt = re.match('chr(.+)', peaknm)
            peak_nm = mt.group(1)
        else:
            peak_nm = peaknm

        final_line = peak_nm + '\t' + col[1] + '\t' + col[2]

        if cellnm in store_target_cell_dic:
            sparse_flt_fl.write(final_line + '\n')

    sparse_fl.close()
    sparse_flt_fl.close()

##updating 081422 not use currently
##we will generate a peak version
def build_peak_version (input_output_dir):

    sparse_flt_fl = input_output_dir + '/opt_downsample_peak.sparse'

    store_peak_dic = {}
    with open (sparse_flt_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            mt = re.match('(.+)_(.+)_(.+)',col[0])
            chrnm = mt.group(1)
            st = mt.group(2)
            ed = mt.group(3)

            final_line = chrnm + '\t' + st + '\t' + ed + '\t' + col[0]
            store_peak_dic[final_line] = 1

    with open (input_output_dir + '/temp_downsample_peak.txt','w+') as opt:
        for eachline in store_peak_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/temp_downsample_peak.txt > ' + \
          input_output_dir + '/temp_downsample_peak_sorted.txt'
    subprocess.call(cmd,shell=True)

##updating 010225
def prepare_motif_acr_data (ipt_R_script,input_motif_acr_fl,input_motif_flt_cutoff,input_peak_tn5_sparse_fl,prefix,input_output_dir):

    cmd = 'Rscript ' + ipt_R_script + \
          ' ' + input_motif_acr_fl + \
          ' ' + input_motif_flt_cutoff + \
          ' ' + input_peak_tn5_sparse_fl + \
          ' ' + prefix + \
          ' ' + input_output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)





















