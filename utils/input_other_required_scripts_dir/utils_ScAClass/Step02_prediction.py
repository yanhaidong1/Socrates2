#!/usr/bin/env python


##this script is to conduct the training and prediction.
##updating 102524 we will check if we will use the low mem


import numpy as np
import sys
import subprocess
import glob
import re
import os
import random


input_training_testing_script_fl = sys.argv[1]

input_training_meta_fl = sys.argv[2]

input_embedding_dir = sys.argv[3]

##parameters
SVM_para = sys.argv[4]

open_smote_upsampling = sys.argv[5]

target_method_string = sys.argv[6]

##save to output dir
input_output_dir = sys.argv[7]

##updating 102524
decide_use_low_mem = sys.argv[8]

##updating 013026
open_parameter_tunning_final = sys.argv[9]



def training_prediction (input_training_testing_script_fl,input_training_meta_fl,input_embedding_dir,
                         SVM_para,open_smote_upsampling,target_method_string,open_parameter_tunning_final,
                         input_output_dir):

    ipt_feature_cell_embeddings_dir = input_embedding_dir

    ipt_PC_dir_list = glob.glob(ipt_feature_cell_embeddings_dir + '/*')

    for eachdir in ipt_PC_dir_list:

        mt = re.match('.+/(.+)',eachdir)
        PCdirnm = mt.group(1)

        mt = re.match('opt_PC(.+)_dir',PCdirnm)
        PCnum = mt.group(1)

        ##this PC is for storing the final output dir
        opt_PC_dir = input_output_dir + '/opt_PC' + PCnum + '_dir'
        if not os.path.exists(opt_PC_dir):
            os.makedirs(opt_PC_dir)

        ipt_training_embeddings_fl = eachdir + '/opt_training_embeddings.csv'
        ipt_testing_embeddings_fl = eachdir + '/opt_indep_testing_embeddings.csv'
        ipt_all_training_meta_fl = input_training_meta_fl

        print('ipt_trainin_meta_fl is')
        print(ipt_all_training_meta_fl)
        print('####')

        print('ipt_training_embeddings_fl is')
        print(ipt_training_embeddings_fl)
        print('####')

        print('ipt_testing_embeddings_fl is')
        print(ipt_testing_embeddings_fl)
        print('####')

        print ('build the one versus all models')

        ##generate a dir in the training dir
        store_one_versus_all_dir = opt_PC_dir + '/store_one_versus_all_dir'
        if not os.path.exists(store_one_versus_all_dir):
            os.makedirs(store_one_versus_all_dir)

        store_split_celltype_dir = opt_PC_dir + '/store_split_celltype_dir'
        if not os.path.exists(store_split_celltype_dir):
            os.makedirs(store_split_celltype_dir)

        store_all_cell_types_dic = {}
        count = 0
        with open (ipt_all_training_meta_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')
                count += 1
                if count != 1:
                    celltype = col[1]
                    store_all_cell_types_dic[celltype] = 1

        for eachcelltype in store_all_cell_types_dic:

            print('the target analyzed celltype is ')
            print(eachcelltype)
            print('####')

            store_eachcelltype_PC_dir = store_one_versus_all_dir + '/store_' + eachcelltype + '_dir'
            if not os.path.exists(store_eachcelltype_PC_dir):
                os.makedirs(store_eachcelltype_PC_dir)

            ##generate new training dataset
            store_final_line_list = []
            first_line = ',cell_type,prob'
            store_final_line_list.append(first_line)
            #count = 0
            with open(ipt_all_training_meta_fl, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')
                    #count += 1
                    #if count != 1:
                    celltype = col[1]

                    if eachcelltype == celltype:

                        #print(celltype)
                        final_line = col[0] + ',' + col[1] + ',' + '1'
                        store_final_line_list.append(final_line)
                    else:
                        newcelltype = 'Others'
                        final_line = col[0] + ',' + newcelltype + ',' + '1'
                        store_final_line_list.append(final_line)
                    #else:
                    #    store_final_line_list.append(eachline)

            with open (store_split_celltype_dir + '/opt_meta_train_' + eachcelltype + '.csv','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            ##now we will do the training
            cmd = 'python ' + input_training_testing_script_fl + \
                  ' ' + ipt_training_embeddings_fl + \
                  ' ' + ipt_testing_embeddings_fl + \
                  ' ' + store_split_celltype_dir + '/opt_meta_train_' + eachcelltype + '.csv' + \
                  ' ' + store_eachcelltype_PC_dir + \
                  ' ' + SVM_para + \
                  ' ' + open_smote_upsampling + \
                  ' ' + target_method_string + \
                  ' ' + open_parameter_tunning_final
            print(cmd)
            subprocess.call(cmd, shell=True)


def training_prediction_low_mem (input_training_testing_script_fl,input_training_meta_fl,input_embedding_dir,
                         SVM_para,open_smote_upsampling,target_method_string,open_parameter_tunning_final,
                         input_output_dir):

    ipt_feature_cell_embeddings_dir = input_embedding_dir

    ipt_PC_dir_list = glob.glob(ipt_feature_cell_embeddings_dir + '/*')

    for eachdir in ipt_PC_dir_list:

        mt = re.match('.+/(.+)',eachdir)
        PCdirnm = mt.group(1)

        mt = re.match('opt_PC(.+)_dir',PCdirnm)
        PCnum = mt.group(1)

        ##this PC is for storing the final output dir
        opt_PC_dir = input_output_dir + '/opt_PC' + PCnum + '_dir'
        if not os.path.exists(opt_PC_dir):
            os.makedirs(opt_PC_dir)

        ##updating 102524
        ##the each dir folder contains the range dir
        ##for each range dir containing the embedding files
        ipt_all_cell_range_dir_list = glob.glob(eachdir + '/*')

        for eachcellrangedir in ipt_all_cell_range_dir_list:

            mt = re.match('.+/(.+)',eachcellrangedir)
            cellrangedirnm = mt.group(1)

            opt_range_dir = opt_PC_dir + '/' + cellrangedirnm
            if not os.path.exists(opt_range_dir):
                os.makedirs(opt_range_dir)

            ipt_training_embeddings_fl = eachcellrangedir + '/opt_training_embeddings.csv'
            ipt_testing_embeddings_fl = eachcellrangedir + '/opt_indep_testing_embeddings.csv'
            ipt_all_training_meta_fl = input_training_meta_fl

            print('ipt_trainin_meta_fl is')
            print(ipt_all_training_meta_fl)
            print('####')

            print('ipt_training_embeddings_fl is')
            print(ipt_training_embeddings_fl)
            print('####')

            print('ipt_testing_embeddings_fl is')
            print(ipt_testing_embeddings_fl)
            print('####')

            print ('build the one versus all models')

            ##generate a dir in the training dir
            store_one_versus_all_dir = opt_range_dir + '/store_one_versus_all_dir'
            if not os.path.exists(store_one_versus_all_dir):
                os.makedirs(store_one_versus_all_dir)

            store_split_celltype_dir = opt_range_dir + '/store_split_celltype_dir'
            if not os.path.exists(store_split_celltype_dir):
                os.makedirs(store_split_celltype_dir)

            store_all_cell_types_dic = {}
            count = 0
            with open (ipt_all_training_meta_fl,'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')
                    count += 1
                    if count != 1:
                        celltype = col[1]
                        store_all_cell_types_dic[celltype] = 1

            for eachcelltype in store_all_cell_types_dic:

                print('the target analyzed celltype is ')
                print(eachcelltype)
                print('####')

                store_eachcelltype_PC_dir = store_one_versus_all_dir + '/store_' + eachcelltype + '_dir'
                if not os.path.exists(store_eachcelltype_PC_dir):
                    os.makedirs(store_eachcelltype_PC_dir)

                ##generate new training dataset
                store_final_line_list = []
                first_line = ',cell_type,prob'
                store_final_line_list.append(first_line)
                #count = 0
                with open(ipt_all_training_meta_fl, 'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split('\t')
                        #count += 1
                        #if count != 1:
                        celltype = col[1]

                        if eachcelltype == celltype:

                            #print(celltype)
                            final_line = col[0] + ',' + col[1] + ',' + '1'
                            store_final_line_list.append(final_line)
                        else:
                            newcelltype = 'Others'
                            final_line = col[0] + ',' + newcelltype + ',' + '1'
                            store_final_line_list.append(final_line)
                        #else:
                        #    store_final_line_list.append(eachline)

                with open (store_split_celltype_dir + '/opt_meta_train_' + eachcelltype + '.csv','w+') as opt:
                    for eachline in store_final_line_list:
                        opt.write(eachline + '\n')

                ##now we will do the training
                cmd = 'python ' + input_training_testing_script_fl + \
                      ' ' + ipt_training_embeddings_fl + \
                      ' ' + ipt_testing_embeddings_fl + \
                      ' ' + store_split_celltype_dir + '/opt_meta_train_' + eachcelltype + '.csv' + \
                      ' ' + store_eachcelltype_PC_dir + \
                      ' ' + SVM_para + \
                      ' ' + open_smote_upsampling + \
                      ' ' + target_method_string + \
                      ' ' + open_parameter_tunning_final
                print(cmd)
                subprocess.call(cmd, shell=True)




if decide_use_low_mem == 'no':

    training_prediction (input_training_testing_script_fl,input_training_meta_fl,input_embedding_dir,
                             SVM_para,open_smote_upsampling,target_method_string,open_parameter_tunning_final,
                             input_output_dir)

else:

    training_prediction_low_mem(input_training_testing_script_fl, input_training_meta_fl, input_embedding_dir,
                                SVM_para, open_smote_upsampling, target_method_string,open_parameter_tunning_final,
                                input_output_dir)





















