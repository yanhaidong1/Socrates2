#!/usr/bin/env python

##this script is to prepare the embeddings from both training and unknown matrix

import sys
import subprocess
import os


input_training_matrix_rds_fl = sys.argv[1]

input_unknown_matrix_csv_fl = sys.argv[2]

input_feature_embedding_R_script_fl = sys.argv[3]

input_training_meta_fl = sys.argv[4]

##other parameter settings
PCnum_string = sys.argv[5]

reduction_major_type = sys.argv[6]

input_output_dir = sys.argv[7]

top_variant_ft_num = sys.argv[8]

pc_variance_cutoff = sys.argv[9]


def prepare_feature_cell_embeddlings (input_training_matrix_rds_fl,input_unknown_matrix_csv_fl,input_training_meta_fl,input_feature_embedding_R_script_fl,
                                      input_output_dir,
                                      PCnum_string,reduction_major_type,top_variant_ft_num,pc_variance_cutoff):


    ##now we only consider the training and independent test dataset
    ##udpating 032722
    PC_list = PCnum_string.split(',')

    for eachPCnum in PC_list:

        opt_PC_dir = input_output_dir + '/opt_PC' + eachPCnum + '_dir'
        if not os.path.exists(opt_PC_dir):
            os.makedirs(opt_PC_dir)

        cmd = 'Rscript ' + input_feature_embedding_R_script_fl + \
              ' ' + input_training_matrix_rds_fl + \
              ' ' + input_unknown_matrix_csv_fl + \
              ' ' + input_training_meta_fl + \
              ' ' + opt_PC_dir + \
              ' ' + eachPCnum + \
              ' ' + reduction_major_type + \
              ' ' + top_variant_ft_num + \
              ' ' + pc_variance_cutoff
        print(cmd)
        subprocess.call(cmd, shell=True)



prepare_feature_cell_embeddlings (input_training_matrix_rds_fl,input_unknown_matrix_csv_fl,input_training_meta_fl,input_feature_embedding_R_script_fl,
                                      input_output_dir,
                                      PCnum_string,reduction_major_type,top_variant_ft_num,pc_variance_cutoff)

















