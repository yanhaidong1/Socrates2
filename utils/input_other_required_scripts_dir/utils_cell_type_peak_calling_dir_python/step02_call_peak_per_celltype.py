#!/usr/bin/env python

import re
import glob
import subprocess
import os
from multiprocessing import Pool
import numpy as np

##updating 062325 fix issue of filtering big peak


def multi_run_step02(args):
    return call_peaks_MACS2_prepare_tn5_peak(*args)

def running_macs2 (ipt_bed_fl,g_size,opt_dir,clustid,pval_or_qval,pqval_num,
                   SLOCAL, LLOCAL, max_gap,
                   extsize,shiftsize):

    ##ori is --extsize 150
    ##shift is -75
    if pval_or_qval == 'qval':

        cmd = 'macs2 callpeak' + \
              ' --cutoff-analysis' + \
              ' -t ' + ipt_bed_fl + \
              ' -f BED ' + \
              ' -g ' + g_size + \
              ' --nomodel ' + \
              ' --keep-dup all ' + \
              ' --extsize ' + extsize + \
              ' --shift ' + shiftsize + \
              ' --max-gap ' + max_gap + \
              ' --slocal ' + SLOCAL + \
              ' --llocal ' + LLOCAL + \
              ' -q ' + pqval_num + \
              ' --outdir ' + opt_dir + \
              ' --bdg ' + \
              ' -n cluster.' + clustid + '.macs2'
        subprocess.call(cmd, shell=True)

    if pval_or_qval == 'pval':
        cmd = 'macs2 callpeak' + \
              ' --cutoff-analysis' + \
              ' -t ' + ipt_bed_fl + \
              ' -f BED ' + \
              ' -g ' + g_size + \
              ' --nomodel ' + \
              ' --keep-dup all ' + \
              ' --extsize ' + extsize + \
              ' --shift ' + shiftsize + \
              ' --max-gap ' + max_gap + \
              ' --slocal ' + SLOCAL + \
              ' --llocal ' + LLOCAL + \
              ' -p ' + pqval_num + \
              ' --outdir ' + opt_dir + \
              ' --bdg ' + \
              ' -n cluster.' + clustid + '.macs2'
        subprocess.call(cmd, shell=True)


##updating 062325
def filter_out_big_peak (ipt_peak_fl,ipt_opt_dir,organ_ct_nm):

    store_final_line_list = []
    with open (ipt_peak_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if (int(col[2]) - int(col[1])) < 1000:
                store_final_line_list.append(eachline)

    with open (ipt_opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def intersect_predicted_peak_bed_fl (peak_fl,final_bed_fl,opt_dir,store_raw_peak_intersect_dir):

    mt = re.match('.+/(.+)', peak_fl)
    flnm = mt.group(1)

    cmd = 'sort -k1,1V -k2,2n ' + peak_fl + ' > ' + peak_fl + '.sorted'
    subprocess.call(cmd, shell=True)

    cmd = 'bedtools intersect -a ' + peak_fl + '.sorted' + ' -b ' + final_bed_fl + ' -c > ' + \
          opt_dir + '/' + flnm + '.sorted.tn5.intersect.txt'
    subprocess.call(cmd, shell=True)

    cmd = 'cp ' + opt_dir + '/' + flnm + '.sorted.tn5.intersect.txt ' +  store_raw_peak_intersect_dir
    subprocess.call(cmd,shell=True)


def permute_random_regions (final_bed_fl,input_gff_fl,peak_fl,opt_dir,input_ref_fai_fl,clustid,store_control_peak_dir,
                            input_unmappable_fl):

    store_excl_bed_line_list = []
    with open (input_gff_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            if not eachline.startswith('#'):

                cate = col[2]
                #if cate == 'gene':
                #    final_line = col[0] + '\t' + col[3] + '\t' + col[4]
                #    store_excl_bed_line_list.append(final_line)

                ##updating 010426 we will also check if it is equal to CDS
                if cate == 'exon' or cate == 'CDS':
                    final_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_excl_bed_line_list.append(final_line)

    with open (peak_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            final_line = col[0] + '\t' + col[1] + '\t' + col[2]
            store_excl_bed_line_list.append(final_line)


    ##updating 040822
    ##exclude the unmappable regions
    if input_unmappable_fl != 'None':
        with open (input_unmappable_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                store_excl_bed_line_list.append(eachline)


    with open (opt_dir + '/opt_' + clustid + '_excl_fl.txt','w+') as opt:
        for eachline in store_excl_bed_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_' + clustid + '_excl_fl.txt > ' + opt_dir + '/opt_' + clustid + '_excl_fl_sorted.txt'
    subprocess.call(cmd,shell=True)

    ##do the shuffle
    cmd = 'bedtools shuffle -i ' + peak_fl + \
          ' -excl ' + opt_dir + '/opt_' + clustid + '_excl_fl_sorted.txt' + \
          ' -g ' + input_ref_fai_fl + \
          ' > ' + opt_dir + '/opt_' + clustid+ '_permute_control_peak.bed'
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_' + clustid + '_permute_control_peak.bed' + \
          ' > ' + opt_dir + '/opt_' + clustid+ '_permute_control_peak_sorted.bed'
    subprocess.call(cmd,shell=True)

    ##do the intersect
    cmd = 'bedtools intersect -a ' + opt_dir + '/opt_' + clustid+ '_permute_control_peak_sorted.bed' + \
          ' -b ' +  final_bed_fl + \
          ' -c > ' + opt_dir + '/opt_' + clustid + '_permute_control_peak_sorted_intersect_tn5.txt'
    subprocess.call(cmd,shell=True)

    cmd = 'cp ' +  opt_dir + '/opt_' + clustid + '_permute_control_peak_sorted_intersect_tn5.txt ' +  store_control_peak_dir
    subprocess.call(cmd,shell=True)


def call_peaks_MACS2_prepare_tn5_peak (input_dir_list,
                      opt_dir,g_size,pval_or_qval,pqval_num,input_ref_fai_fl,input_gff_fl,
                                       store_MACS2_raw_peak_dir,store_control_peak_dir,store_raw_peak_intersect_dir,
                                       SLOCAL, LLOCAL, max_gap,extsize,shiftsize,input_unmappable_fl,s2_open_callpeak,s2_open_permute):

    for eachdir in input_dir_list:

        ##check whether the dir is file
        ##if it is dir, it means we will use the bed file of rice version
        if os.path.isfile(eachdir) != True:

            #print(eachdir)

            mt = re.match('.+/(.+)',eachdir)
            dirnm = mt.group(1)

            ##extract the name for the cell type name
            ##becaureful about the pattern
            ##it should be crownroot.crownroot.Atrichoblast

            ##updating 072923 not organ mode only works for the no organ mode
            #if s2_open_rice_organ_mode == 'yes':
            mt = re.match('.+\.(.+\..+)',dirnm)
            organ_ct_nm = mt.group(1)
            ##crownroot.Atrichoblast
            #else:


            fl_list = glob.glob(eachdir + '/*')

            ##check whether pool is in the flnm
            pool_fl_num = 0
            for eachfl in fl_list:
                mt = re.match('.+/(.+)', eachfl)
                flnm = mt.group(1)

                if 'pool' in flnm:
                    pool_fl_num += 1


            for eachfl in fl_list:

                #print(eachfl)

                mt = re.match('.+/(.+)',eachfl)
                flnm = mt.group(1)


                if pool_fl_num != 0:

                    if 'pool' in flnm:

                        final_bed_fl = eachfl


                        if s2_open_callpeak == 'yes':
                            ##do the calling peaks
                            running_macs2(final_bed_fl, g_size, opt_dir, organ_ct_nm, pval_or_qval, pqval_num,SLOCAL, LLOCAL, max_gap,extsize,shiftsize)

                            ##intersect with predicted peaks
                            temp_peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.narrowPeak'
                            filter_out_big_peak(temp_peak_fl, opt_dir, organ_ct_nm)
                            peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'

                            cmd = 'cp ' + peak_fl + ' ' + store_MACS2_raw_peak_dir
                            subprocess.call(cmd,shell=True)
                            intersect_predicted_peak_bed_fl(peak_fl, final_bed_fl, opt_dir,store_raw_peak_intersect_dir)

                        if s2_open_permute == 'yes':
                            ##updating 032422 we will continue for the function
                            ##shuffle
                            peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'
                            permute_random_regions(final_bed_fl, input_gff_fl, peak_fl, opt_dir, input_ref_fai_fl, organ_ct_nm,store_control_peak_dir,input_unmappable_fl)

                else:
                    final_bed_fl = eachfl
                    #print(final_bed_fl)

                    if s2_open_callpeak == 'yes':
                        ##do the calling peaks
                        running_macs2(final_bed_fl, g_size, opt_dir, organ_ct_nm, pval_or_qval, pqval_num, SLOCAL, LLOCAL,
                                      max_gap,extsize,shiftsize)

                        ##intersect with predicted peaks
                        temp_peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.narrowPeak'
                        filter_out_big_peak(temp_peak_fl, opt_dir, organ_ct_nm)
                        peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'

                        cmd = 'cp ' + peak_fl + ' ' + store_MACS2_raw_peak_dir
                        subprocess.call(cmd, shell=True)
                        intersect_predicted_peak_bed_fl(peak_fl, final_bed_fl, opt_dir, store_raw_peak_intersect_dir)

                    ##updating 032422 we will continue for the function
                    ##shuffle
                    if s2_open_permute == 'yes':
                        peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'

                        permute_random_regions(final_bed_fl, input_gff_fl, peak_fl, opt_dir, input_ref_fai_fl, organ_ct_nm,
                                               store_control_peak_dir, input_unmappable_fl)

        else:
            ##if the file is bed file we will think it is from the soybean or other format
            final_bed_fl = eachdir

            mt = re.match('.+/(.+)',final_bed_fl)
            flnm = mt.group(1)
            #mt = re.match('cluster_(.+)\.Gm_seed_atlas.bed',flnm)
            ##updating 101422
            mt = re.match('opt_(.+)\.bed',flnm)
            organ_ct_nm = mt.group(1)

            if s2_open_callpeak == 'yes':
                ##do the calling peaks
                running_macs2(final_bed_fl, g_size, opt_dir, organ_ct_nm, pval_or_qval, pqval_num, SLOCAL, LLOCAL,
                              max_gap,extsize,shiftsize)

                ##intersect with predicted peaks
                temp_peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.narrowPeak'
                filter_out_big_peak(temp_peak_fl, opt_dir, organ_ct_nm)
                peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'

                cmd = 'cp ' + peak_fl + ' ' + store_MACS2_raw_peak_dir
                print(cmd)
                subprocess.call(cmd, shell=True)
                intersect_predicted_peak_bed_fl(peak_fl, final_bed_fl, opt_dir, store_raw_peak_intersect_dir)

            ##updating 032422 we will continue for the function
            ##shuffle
            if s2_open_permute == 'yes':

                peak_fl = opt_dir + '/cluster.' + organ_ct_nm + '.macs2_peaks.flt.narrowPeak'

                permute_random_regions(final_bed_fl, input_gff_fl, peak_fl, opt_dir, input_ref_fai_fl, organ_ct_nm,
                                       store_control_peak_dir, input_unmappable_fl)



def run_pipeline (input_output_dir,ipt_target_bed_fl_list,input_gff_fl,input_ref_fai_fl,
                  input_unmappable_fl,
                  input_core_num,
                  g_size,
                  pval_or_qval,
                  pqval_num,
                  SLOCAL, LLOCAL, max_gap,extsize,shiftsize,
                  open_callpeak,open_permute):


    store_MACS2_calling_dir = input_output_dir + '/store_MACS2_calling_dir'
    if not os.path.exists(store_MACS2_calling_dir):
        os.makedirs(store_MACS2_calling_dir)

    store_MACS2_raw_peak_dir = input_output_dir + '/store_MACS2_raw_peak_dir'
    if not os.path.exists(store_MACS2_raw_peak_dir):
        os.makedirs(store_MACS2_raw_peak_dir)

    store_control_peak_dir = input_output_dir + '/store_control_peak_intersect_dir'
    if not os.path.exists(store_control_peak_dir):
        os.makedirs(store_control_peak_dir)

    store_raw_peak_intersect_dir = input_output_dir + '/store_raw_peak_intersect_dir'
    if not os.path.exists(store_raw_peak_intersect_dir):
        os.makedirs(store_raw_peak_intersect_dir)

    ##now it does not work for the header
    ##store the cluster id
    #if header == '':
    #    pool_bed_dir_list = glob.glob(input_pool_bed_dir + '/*')
    #else:
    #    pool_bed_dir_list = glob.glob(input_pool_bed_dir + '/' + header + '*')

    pool_bed_dir_list = ipt_target_bed_fl_list
    print('analyzed pool bed dir list is')
    print(pool_bed_dir_list)

    #print('pool bed dir list is')
    #print(pool_bed_dir_list)

    ipt_number_bed_dir = len(pool_bed_dir_list)

    if ipt_number_bed_dir <= int(input_core_num):
        input_core_num_final = ipt_number_bed_dir
    else:

        input_core_num_final = int(input_core_num)


    array_split = np.array_split(pool_bed_dir_list, int(input_core_num_final))




    for x in range(0, int(input_core_num_final)):
        dir_code = x + 1
        temp_output_dir = store_MACS2_calling_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)
    temp_all_output_dir_list = glob.glob(store_MACS2_calling_dir + '/*')

    ##store the core target bed file list
    store_core_dic = {}
    core_count = -1
    for eacharray in array_split:
        core_count += 1

        ##the array list contains cluster information
        array_list = []
        for eachitem in eacharray:
            array_list.append(eachitem)

        store_core_dic[str(core_count)] = array_list


    ##generate a list to store the define_hetero function
    run_list = []
    for x in range(0, int(input_core_num_final)):

        print(store_core_dic[str(x)])
        print(temp_all_output_dir_list[x])

        each_func_argument = (store_core_dic[str(x)],
                              temp_all_output_dir_list[x],
                              g_size,pval_or_qval,pqval_num,
                              input_ref_fai_fl,input_gff_fl,
                              store_MACS2_raw_peak_dir,store_control_peak_dir,store_raw_peak_intersect_dir,
                              SLOCAL, LLOCAL, max_gap,extsize,shiftsize,
                              input_unmappable_fl,open_callpeak,open_permute
                              )
        run_list.append(each_func_argument)

    pool = Pool(int(input_core_num_final))
    pool.map(multi_run_step02, run_list)


