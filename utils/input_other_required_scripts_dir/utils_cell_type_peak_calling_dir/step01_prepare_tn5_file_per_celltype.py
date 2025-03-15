#!/usr/bin/env python

import glob
import subprocess
import os
from multiprocessing import Pool

def multi_run_step01(args):
    return prepare_bed(*args)

def prepare_bed (input_tn5_fl,input_target_cluster_fl,input_target_meta_fl,target_colnm,opt_dir):

    store_target_clusternm_list = []
    with open (input_target_cluster_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            store_target_clusternm_list.append(eachline)

    for eachclusternm in store_target_clusternm_list:

        store_target_cellid_dic = {}
        count = 0
        target_index = 0
        with open (input_target_meta_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1

                if count == 1:

                    #target_index = col.index(target_colnm)

                    temp_target_index = col.index(target_colnm)

                    all_len_col = len(col)

                    target_index = temp_target_index - all_len_col


                else:

                    # target_cluster_col = col[int(s1_targetclust_colNum) - 1]
                    # target_cluster_col = col[-1]
                    target_cluster_col = col[int(target_index)]
                    if target_cluster_col == eachclusternm:
                        store_target_cellid_dic[col[0]] = 1



        bed_file = open(input_tn5_fl,'r')
        target_bed_file = open (opt_dir + '/opt_' + eachclusternm + '.bed','w')

        for eachline in bed_file:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cellnm = col[3]

            if cellnm in store_target_cellid_dic:
                target_bed_file.write(eachline + '\n')

        bed_file.close()
        target_bed_file.close()


def prepare_bed_files_parallele (input_tn5_fl, input_meta_fl,target_colnm,
                           input_core_num,
                           input_output_dir):


    print ('- step01 build bed files across all clusters')

    step01_prepare_bed_file_dir = input_output_dir

    store_pool_bed_dir = step01_prepare_bed_file_dir + '/store_pool_bed_dir'
    if not os.path.exists(store_pool_bed_dir):
        os.makedirs(store_pool_bed_dir)

    store_target_meta_fl_dir = step01_prepare_bed_file_dir + '/store_target_meta_fl_dir'
    if not os.path.exists(store_target_meta_fl_dir):
        os.makedirs(store_target_meta_fl_dir)

    store_temp_cluster_fl_dir = step01_prepare_bed_file_dir + '/store_temp_cluster_fl_dir'
    if not os.path.exists(store_temp_cluster_fl_dir):
        os.makedirs(store_temp_cluster_fl_dir)


    ########
    ##step02 using the target file, we will build the different dirs
    target_cluster_dic = {}
    count = 0

    ##updating 031525
    target_index = 0

    with open (input_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            #target_cluster_col = col[int(s1_targetclust_colNum)-1]
            #target_cluster_col = col[-1]
            count += 1
            if count == 1:

                temp_target_index = col.index(target_colnm)

                all_len_col = len(col)

                target_index = temp_target_index - all_len_col

            else:

                target_cluster_col = col[int(target_index)]

                target_cluster_dic[target_cluster_col] = 1






    target_cluster_list = list(target_cluster_dic.keys())

    with open (store_target_meta_fl_dir + '/temp_target_cluster_list_fl.txt','w+') as opt:
        for eachline in target_cluster_list:
            opt.write(eachline + '\n')

    temp_all_split_cluster_dir = step01_prepare_bed_file_dir + '/temp_all_split_cluster_dir'
    if not os.path.exists(temp_all_split_cluster_dir):
        os.makedirs(temp_all_split_cluster_dir)

    ##create target chr fl
    cmd = 'split -n l/' + str(input_core_num) + ' ' + store_target_meta_fl_dir + '/temp_target_cluster_list_fl.txt' + ' ' + temp_all_split_cluster_dir + '/temp_split_'
    subprocess.call(cmd,shell=True)


    ##create output dir
    for x in range(0, int(input_core_num)):
        dir_code = x + 1
        temp_output_dir = store_pool_bed_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    temp_split_file_list = glob.glob(temp_all_split_cluster_dir + '/temp_split*')
    temp_all_output_dir_list = glob.glob(store_pool_bed_dir + '/*')


    ########
    ##step03 run the pipeline
    ##prepare_bed (input_subset_tn5_fl,input_target_cluster_fl,input_target_meta_fl,opt_dir,
    ##             s1_targetclust_colNum)
    pool = Pool(int(input_core_num))
    run_list = []
    for x in range(0, int(input_core_num)):
        each_func_argument = (input_tn5_fl,
                              temp_split_file_list[x],
                              input_meta_fl,target_colnm,
                              temp_all_output_dir_list[x])
        run_list.append(each_func_argument)
    pool.map(multi_run_step01, run_list)

