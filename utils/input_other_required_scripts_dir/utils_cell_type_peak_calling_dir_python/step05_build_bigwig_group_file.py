#!/usr/bin/env python

##this script we will build the bigwig for parallel

#import argparse
import subprocess
import re
import glob
#import sys
import os
from multiprocessing import Pool
import numpy as np



def multi_run_step05(args):
    return build_bigwig(*args)

def check_read_num (ipt_bed_fl):
    count = 0
    with open (ipt_bed_fl,'r') as ipt:
        for eachline in ipt:
            count += 1
    return (count)

def build_bigwig (input_required_script_dir,ipt_meta_fl,ipt_bed_fl_list,input_ref_fl,
                  ipt_opt_dir,opt_final_dir,
                  g_size,extsize,shiftsize,max_gap,SLOCAL,LLOCAL):

    #print(ipt_bed_fl_list)

    for eachbed_fl in ipt_bed_fl_list:

        mt = re.match('.+/(.+)',eachbed_fl)
        flnm = mt.group(1)
        mt = re.match('opt_(.+)\.bed',flnm)
        celltype = mt.group(1)

        store_celltype_ctgroup_dic = {}
        count = 0
        with open (ipt_meta_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    group = col[-1]
                    celltype_meta = col[-3]
                    cellID = col[0]

                    if celltype_meta in store_celltype_ctgroup_dic:
                        if group in store_celltype_ctgroup_dic[celltype_meta]:
                            store_celltype_ctgroup_dic[celltype_meta][group].append(cellID)
                        else:
                            store_celltype_ctgroup_dic[celltype_meta][group]= []
                            store_celltype_ctgroup_dic[celltype_meta][group].append(cellID)
                    else:
                        store_celltype_ctgroup_dic[celltype_meta] = {}
                        store_celltype_ctgroup_dic[celltype_meta][group] = []
                        store_celltype_ctgroup_dic[celltype_meta][group].append(cellID)

        #print(store_celltype_ctgroup_dic)

        for eachgroup in store_celltype_ctgroup_dic[celltype]:

            target_cell_id_list = store_celltype_ctgroup_dic[celltype][eachgroup]

            celltype_group_id = celltype + '__' + eachgroup

            ##preform the filtration


            bed_file = open(eachbed_fl, 'r')
            target_bed_file = open(ipt_opt_dir + '/opt_' + celltype_group_id + '.bed', 'w')

            for eachline in bed_file:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')
                cellnm = col[3]

                if cellnm in target_cell_id_list:
                    target_bed_file.write(eachline + '\n')

            bed_file.close()
            target_bed_file.close()

            pqval_num = '0.05'


            ##perform the peak calling
            cmd = 'macs2 callpeak' + \
                  ' --cutoff-analysis' + \
                  ' -t ' + ipt_opt_dir + '/opt_' + celltype_group_id + '.bed' + \
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
                  ' --outdir ' + ipt_opt_dir + \
                  ' --bdg ' + \
                  ' -n cluster.' + celltype_group_id + '.macs2'
            print(cmd)
            subprocess.call(cmd, shell=True)


            target_bdg_fl = ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.bdg'
            cleanBED_script = input_required_script_dir + '/cleanBED.py'
            readdepth = check_read_num(ipt_opt_dir + '/opt_' + celltype_group_id + '.bed')

            if os.path.exists(ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg'):

                print('the file ' + ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg' + ' exists')

                ##check file is empty or not

                if os.stat(ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg').st_size != 0:
                    print(
                        'the file ' + ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg' + ' is not empty')

                else:
                    ##generate bedgrapth
                    cmd = 'sort -k1,1 -k2,2n ' + target_bdg_fl + ' | ' + \
                          ' python ' + cleanBED_script + ' ' + input_ref_fl + ' ' + str(readdepth) + ' - >' + \
                          ' ' + ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg'
                    print(cmd)
                    subprocess.call(cmd, shell=True)

            else:

                ##generate bedgrapth
                cmd = 'sort -k1,1 -k2,2n ' + target_bdg_fl + ' | ' + \
                      ' python ' + cleanBED_script + ' ' + input_ref_fl + ' ' + str(readdepth) + ' - >' + \
                      ' ' + ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg'
                print(cmd)
                subprocess.call(cmd, shell=True)

            cmd = 'bedGraphToBigWig ' + ipt_opt_dir + '/cluster.' + celltype_group_id + '.macs2_treat_pileup.clean.bdg ' + input_ref_fl + ' ' + opt_final_dir + '/' + celltype_group_id + '.bw'
            print(cmd)
            subprocess.call(cmd, shell=True)





def build_bigwig_fl_parallele (input_output_dir,
                               step05_build_bigwig_fl_dir,input_required_script_dir,
                               input_ref_fl,ipt_meta_fl,input_core_num,
                               g_size, extsize, shiftsize, max_gap, SLOCAL, LLOCAL):

    store_final_bw_dir = step05_build_bigwig_fl_dir + '/store_final_bw_dir'
    if not os.path.exists(store_final_bw_dir):
        os.makedirs(store_final_bw_dir)


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
    #store_store_final_bw_dir = step04_build_bigwig_fl_dir + '/store_final_bw_dir'
    #if not os.path.exists(store_store_final_bw_dir):
    #    os.makedirs(store_store_final_bw_dir)

    #store_store_final_bdg_dir = step04_build_bigwig_fl_dir + '/store_final_bdg_dir'
    #if not os.path.exists(store_store_final_bdg_dir):
    #    os.makedirs(store_store_final_bdg_dir)


    store_temp_opt_fl_dir = step05_build_bigwig_fl_dir + '/store_temp_opt_fl_dir'
    if not os.path.exists(store_temp_opt_fl_dir):
        os.makedirs(store_temp_opt_fl_dir)

    ##updating 031525 build the exact number of file
    store_alltargetOrgan_bed_fl_num = len(store_alltargetOrgan_bed_fl_list)

    if store_alltargetOrgan_bed_fl_num <= int(input_core_num):
        input_core_num_final = store_alltargetOrgan_bed_fl_num
    else:
        input_core_num_final = int(input_core_num)


    ##split the target bed file into different
    store_core_dic = {}
    core_count = -1
    array_split_list = np.array_split(store_alltargetOrgan_bed_fl_list, int(input_core_num_final))
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
    for x in range(0, int(input_core_num_final)):
        dir_code = x + 1
        temp_output_dir = store_temp_opt_fl_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    temp_all_output_dir_list = glob.glob(store_temp_opt_fl_dir + '/*')

    ##build_bigwig_fl (ipt_bed_fl_list,input_output_dir,input_required_script_dir,
    #                 input_ref_fl,ipt_opt_dir,
    #                 store_final_bw_dir)

    #s2_open_rice_organ_mode = 'no'

    pool = Pool(int(input_core_num_final))
    run_list = []
    for x in range(0, int(input_core_num_final)):
        each_func_argument = (input_required_script_dir,
                              ipt_meta_fl,
                              store_core_dic[x],
                              input_ref_fl,
                              temp_all_output_dir_list[x],
                              store_final_bw_dir,
                              g_size, extsize, shiftsize, max_gap, SLOCAL, LLOCAL)

        run_list.append(each_func_argument)
    pool.map(multi_run_step05, run_list)

    ##build_bigwig (input_required_script_dir,ipt_meta_fl,ipt_bed_fl_list,input_ref_fl,
    #              ipt_opt_dir,opt_final_dir,
    #              g_size,extsize,shiftsize,max_gap,SLOCAL,LLOCAL)

    #build_bigwig(input_required_script_dir, ipt_meta_fl, ipt_bed_fl_list, input_ref_fl,
    #             ipt_opt_dir, opt_final_dir,
    #             g_size, extsize, shiftsize, max_gap, SLOCAL, LLOCAL)