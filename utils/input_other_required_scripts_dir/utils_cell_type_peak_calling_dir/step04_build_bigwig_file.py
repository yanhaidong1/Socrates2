#!/usr/bin/env python

##this step is to build the bigwig file that could be used for the browser investigating.
import argparse
import subprocess
import re
import glob
import sys
import os
from multiprocessing import Pool
import numpy as np



def multi_run_step04(args):
    return build_bigwig_fl(*args)

def check_read_num (ipt_bed_fl):
    count = 0
    with open (ipt_bed_fl,'r') as ipt:
        for eachline in ipt:
            count += 1
    return (count)

##updating 111422
def check_cell_num (ipt_bed_fl):

    store_cellcount_dic = {}
    with open (ipt_bed_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cellnm = col[3]
            store_cellcount_dic[cellnm] = 1
    count = len(list(store_cellcount_dic.keys()))
    return (count)

def build_bigwig_fl (ipt_bed_fl_list,input_output_dir,input_required_script_dir,
                     input_ref_fl,ipt_opt_dir,
                     store_final_bw_dir,
                     store_store_final_bdg_dir,s2_open_rice_organ_mode,
                     s4_open_sort_bdg,
                     s4_norm_type):

    for eachbedfl in ipt_bed_fl_list:

        mt = re.match('.+/(.+)',eachbedfl)
        flnm = mt.group(1)

        mt = re.match('opt_(.+)\.bed', flnm)
        organct = mt.group(1)

        if s2_open_rice_organ_mode == 'yes':

            mt = re.match('opt_(.+)\..+\.bed', flnm)
            organnm = mt.group(1)

        else:

            organnm = 'all'

        print ('target analyzed organct is ' + organct)

        ##we will find the corresponding fl in the calling peak folder
        #step02_call_peaks_dir = input_output_dir + '/s2_open_call_peak_final_dir'

        #valdir_list = glob.glob(step02_call_peaks_dir + '/' + organnm + '/*')
        #valdir = valdir_list[0]

        alltemp_dir_list = glob.glob(input_output_dir + '/s2_open_call_peak_final_dir' + '/store_MACS2_calling_dir/*')


        store_all_organct_bdg_fl_dic = {}
        for eachtemdir in alltemp_dir_list:
            ## opt_dir + '/cluster.' + clustid + '.macs2_treat_pileup.bdg
            allfl_list = glob.glob(eachtemdir + '/*')
            for eachfl in allfl_list:
                mt = re.match('.+/(.+)',eachfl)
                flnm = mt.group(1)
                if re.match('cluster\.(.+)\.macs2_treat_pileup.bdg',flnm):

                    target_bdg_fl = eachfl
                    mt = re.match('cluster\.(.+)\.macs2_treat_pileup.bdg',flnm)
                    organct_nm = mt.group(1)

                    store_all_organct_bdg_fl_dic[organct_nm] = target_bdg_fl

        print ('store_all_organct_bdg_fl_dic is ')
        print (store_all_organct_bdg_fl_dic)

        target_bdg_fl = store_all_organct_bdg_fl_dic[organct]

        print ('target organct is ' + organct)
        print ('target analyzed bdg fl is')
        print (target_bdg_fl)

        readdepth = ''
        cleanBED_script = ''
        if s4_norm_type == 'tn5count':
            readdepth = check_read_num(eachbedfl)
            cleanBED_script = input_required_script_dir + '/cleanBED.pl'
        else:
            if s4_norm_type == 'cellcount':
                readdepth = check_cell_num(eachbedfl)
                cleanBED_script = input_required_script_dir + '/cleanBED_cellpropWay.pl'

            #else:
                ##udpating 012523
                #if s4_norm_type == 'peaktn5count':

                    ##we should find the peak first
                #    store_organ_pqval_num_dic = 'no'
                #    pqval_cutoff = store_organ_pqval_num_dic[organnm]
                #    fdr_choose = finalfdr_list[-1]
                #    target_peak_fl = input_output_dir + '/step02_call_peaks_dir/' + organnm + '/' + pqval_cutoff + '/' + fdr_choose + '/store_filtered_peak_dir/opt_flt_' + organct + '_peak.txt'

                    ##we will check if the file existing
                #    if os.path.isfile(target_peak_fl) == True:
                #        print('the target peak fl exists')
                #    else:
                #        print('please make sure the target peak fl is existing')
                #        break

                    ##do the intersecting
                #    cmd = 'bedtools intersect -a ' + eachbedfl + ' -b ' + target_peak_fl + ' -wa -u > ' + ipt_opt_dir + '/temp_' + organct + '_peak_tn5.bed'
                #    print(cmd)
                #    subprocess.call(cmd,shell=True)

                #    readdepth = check_read_num(ipt_opt_dir + '/temp_' + organct + '_peak_tn5.bed')
                #    cleanBED_script = input_required_script_dir + '/cleanBED.pl'

                #else:
                #    print('make sure we use cellcount or tn5count other than other words')
                #    break



        ##updating 102022 check whether the organct already existed
        if os.path.exists(ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg'):

            print('the file ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg' + ' exists')

            ##check file is empty or not

            if os.stat(ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg').st_size != 0:
                print('the file ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg' + ' is not empty')

            else:
                ##generate bedgrapth
                cmd = 'sort -k1,1 -k2,2n ' + target_bdg_fl + ' | ' + \
                      ' perl ' + cleanBED_script + ' ' + input_ref_fl + ' ' + str(readdepth) + ' - >' + \
                      ' ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg'
                print(cmd)
                subprocess.call(cmd, shell=True)

        else:
            ##generate bedgrapth
            cmd = 'sort -k1,1 -k2,2n ' + target_bdg_fl + ' | ' + \
                  ' perl ' + cleanBED_script + ' ' + input_ref_fl + ' ' + str(readdepth) + ' - >' + \
                  ' ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg'
            print(cmd)
            subprocess.call(cmd, shell=True)


        cmd = 'bedGraphToBigWig ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg ' + input_ref_fl + ' ' + store_final_bw_dir + '/cluster.' + organct + '.raw.bw'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##updating 101922 check whether we need to sort the bdg file
        if s4_open_sort_bdg == 'yes':
            cmd = 'sort -k1,1V -k2,2n ' + ipt_opt_dir + '/cluster.' + organct + '.macs2_treat_pileup.clean.bdg > ' + store_store_final_bdg_dir + '/cluster_' + organct + '.' + organnm + '.raw.sorted.bdg'
            print(cmd)
            subprocess.call(cmd,shell=True)




def build_bigwig_fl_parallele (input_output_dir,input_final_output_dir,input_required_script_dir,input_ref_fl,input_core_num,
                     s4_open_sort_bdg,s4_norm_type):

    step04_build_bigwig_fl_dir = input_final_output_dir

    input_pool_bed_dir = input_output_dir + '/s1_open_prepare_tn5_file_final_dir/store_pool_bed_dir/'
    ##we first need to build a dic to store all the target organ bed file
    all_temp_output_dir_list = glob.glob(input_pool_bed_dir + '/*')
    store_all_bed_fl_path_list = []
    for eachdir in all_temp_output_dir_list:
        all_bed_fl_list = glob.glob(eachdir + '/*')
        for eachbedfl in all_bed_fl_list:
            store_all_bed_fl_path_list.append(eachbedfl)

    #with open (input_output_dir + '/debug_all_bed_fl_path.txt','w+') as opt:
    #    for eachline in store_all_bed_fl_path_list:
    #        opt.write(eachline + '\n')

    #s4_target_organ_list = s4_target_organ_string.split(',')

    store_alltargetOrgan_bed_fl_list = store_all_bed_fl_path_list

    ##do the parallele
    store_store_final_bw_dir = step04_build_bigwig_fl_dir + '/store_store_final_bw_dir'
    if not os.path.exists(store_store_final_bw_dir):
        os.makedirs(store_store_final_bw_dir)

    store_store_final_bdg_dir = step04_build_bigwig_fl_dir + '/store_store_final_bdg_dir'
    if not os.path.exists(store_store_final_bdg_dir):
        os.makedirs(store_store_final_bdg_dir)


    store_temp_opt_fl_dir = step04_build_bigwig_fl_dir + '/store_temp_opt_fl_dir'
    if not os.path.exists(store_temp_opt_fl_dir):
        os.makedirs(store_temp_opt_fl_dir)


    ##split the target bed file into different
    store_core_dic = {}
    core_count = -1
    array_split_list = np.array_split(store_alltargetOrgan_bed_fl_list, int(input_core_num))
    for eacharray in array_split_list:
        core_count += 1

        ##the array list contains cluster information
        array_list = []
        for eachitem in eacharray:
            array_list.append(eachitem)

        store_core_dic[core_count] = array_list

    print('the input file list is')
    print(store_core_dic)

    ##create output dir
    for x in range(0, int(input_core_num)):
        dir_code = x + 1
        temp_output_dir = store_temp_opt_fl_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    temp_all_output_dir_list = glob.glob(store_temp_opt_fl_dir + '/*')

    ##build_bigwig_fl (ipt_bed_fl_list,input_output_dir,input_required_script_dir,
    #                 input_ref_fl,ipt_opt_dir,
    #                 store_final_bw_dir)

    s2_open_rice_organ_mode = 'no'

    pool = Pool(int(input_core_num))
    run_list = []
    for x in range(0, int(input_core_num)):
        each_func_argument = (store_core_dic[x],
                              input_output_dir,
                              input_required_script_dir,
                              input_ref_fl,
                              temp_all_output_dir_list[x],
                              store_store_final_bw_dir,
                              store_store_final_bdg_dir,
                              s2_open_rice_organ_mode,
                              s4_open_sort_bdg,
                              s4_norm_type)
        run_list.append(each_func_argument)
    pool.map(multi_run_step04, run_list)




