#!/usr/bin/env python

import argparse
import glob
import sys
import os
import subprocess
import re


##this script is to conduct the clustering after the loading the data

def get_parsed_args():

    parser = argparse.ArgumentParser(description="Load data and perform QC")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    parser.add_argument("-ipt_obj_dir", dest='input_object_dir', help = 'Users need to provide a directory that stores the libraries object obtained from the loading data and quality control step.')

    ###############################
    ##optional essential parameters
    parser.add_argument('-num_var', dest = 'number_variable_windows', help = 'Number of variable windows used for the clustering.'
                                                                             'Default: 20000')

    parser.add_argument('-dim_type', dest = 'dimension_reduction_type', help = 'Type of dimension reduction for the data. Two options: NMF or SVD.'
                                                                               'Default: SVD')

    parser.add_argument('-num_PCs', dest = 'number_PCs', help = 'Number of PCs used for clustering.'
                                                                'Default: 30')

    parser.add_argument('-norm', dest = 'norm_method', help = 'Normalization method for the data matrix. Two options: tfidf or regmodel.'
                                                              'Default: tfidf')

    parser.add_argument('-res', dest = 'cluster_resolution', help = 'Resolution of clustering.'
                                                                    'Default: 0.5')

    parser.add_argument('-harmony', dest = 'open_harmony', help = 'Perform harmony for data from different libraries.'
                                                                  'Default: yes')

    parser.add_argument('-rm_temp', dest = 'remove_temp_file', help = 'Remove the temp files built by each step.'
                                                                      'Default: yes')

    #####################
    ##optional parameters
    parser.add_argument('-min_cell', dest = 'minimum_cell_num', help = 'Minimum number of accessible features per cell for cell filtration.'
                                                                       'Default: 1000')

    parser.add_argument('-min_ft_freq', dest = 'minimum_feature_freq', help = 'Minimum peak accessibility frequency lower quantile cut-offs.'
                                                                              'Default: 0.005')

    parser.add_argument('-max_ft_freq', dest = 'maximum_feature_freq', help = 'Maximum peak accessibility frequency upper quantile cut-offs.'
                                                                              'Default: 0.005')

    parser.add_argument('-core_doublet', dest = 'core_num_doublet', help = 'Set a core number for doublet calling.'
                                                                                          'Default: 1')


    parser.add_argument('-only_cluster', dest = 'open_only_cluster' , help = 'Only open the clustering step. Users could adjust the resolution to obtain different clustering performance.'
                                                                             'Default: no')



    ##parse of parameters
    args = parser.parse_args()
    return args



def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    if args.required_script_dir is None:
        print('Cannot find required script dir, please provide the dir in \'-script_dir\' !')
        return

    input_required_scripts_dir = args.required_script_dir

    store_final_parameter_line_list = []

    if args.input_object_dir is None:
        print('Cannot find input object dir, please provide the dir in \'-ipt_obj_dir\' !')
        return

    ##step01 parameters
    if args.number_variable_windows is not None:
        num_var_final = args.number_variable_windows
    else:
        num_var_final = '20000'

    store_final_parameter_line_list.append('num_var_final <- ' + num_var_final)

    if args.dimension_reduction_type is not None:
        rdType = args.dimension_reduction_type
    else:
        rdType = 'SVD'

    store_final_parameter_line_list.append('rdType <- ' + '\'' + rdType + '\'')


    if args.number_PCs is not None:
        numPCs = args.number_PCs
    else:
        numPCs = '30'

    store_final_parameter_line_list.append('numPCs <- ' + numPCs)

    if args.norm_method is not None:
        normalize_way = args.norm_method
    else:
        normalize_way = 'tfidf'

    store_final_parameter_line_list.append('normalize_way <- ' + '\'' + normalize_way + '\'')

    if args.cluster_resolution is not None:
        res_clust = args.cluster_resolution
    else:
        res_clust = '0.5'

    store_final_parameter_line_list.append('res_clust <- ' + res_clust)


    if args.open_harmony is not None:
        do_harmony = args.open_harmony
    else:
        do_harmony = 'yes'

    store_final_parameter_line_list.append('do_harmony <- ' + '\'' + do_harmony + '\'')


    if args.remove_temp_file is not None:
        remove_temp = args.remove_temp_file
    else:
        remove_temp = 'yes'

    store_final_parameter_line_list.append('remove_temp <- ' + '\'' + remove_temp + '\'')


    if args.minimum_cell_num is not None:
        min_cell_final = args.minimum_cell_num
    else:
        min_cell_final = '1000'

    store_final_parameter_line_list.append('min_cell_final <- ' + min_cell_final)

    if args.minimum_feature_freq is not None:
        min_ft_freq_final = args.minimum_feature_freq
    else:
        min_ft_freq_final = '0.005'

    store_final_parameter_line_list.append('min_ft_freq_final <- ' + min_ft_freq_final)

    if args.maximum_feature_freq is not None:
        max_ft_freq_final = args.maximum_feature_freq
    else:
        max_ft_freq_final = '0.005'

    store_final_parameter_line_list.append('max_ft_freq_final <- ' + max_ft_freq_final)

    if args.core_num_doublet is not None:
        core_doublet_final = args.core_num_doublet
    else:
        core_doublet_final = '1'

    store_final_parameter_line_list.append('core_doublet_final <- ' + core_doublet_final)

    if args.open_only_cluster is not None:
        only_cluster = args.open_only_cluster
    else:
        only_cluster = 'no'

    store_final_parameter_line_list.append('only_cluster <- ' + only_cluster)

    ##we will firstly build up the parameter setting file
    with open(output_dir + '/temp_defined_parameters.config', 'w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')

    ##obtain the previous folder
    ipt_path_to_preload_R = ''

    if re.match('(.+)/utils/input_other_required_scripts_dir/', input_required_scripts_dir):
        mt = re.match('(.+)/utils/input_other_required_scripts_dir/', input_required_scripts_dir)
        ipt_path_to_preload_R = mt.group(1) + '/R'

    if re.match('(.+)/utils/input_other_required_scripts_dir', input_required_scripts_dir):
        mt = re.match('(.+)/utils/input_other_required_scripts_dir', input_required_scripts_dir)
        ipt_path_to_preload_R = mt.group(1) + '/R'


    ipt_script = input_required_scripts_dir + '/utils_load_and_QC_data/clustering.R'
    ipt_obj_dir = args.input_object_dir
    ipt_config_fl = output_dir + '/temp_defined_parameters.config'

    ##running
    cmd = 'Rscript ' + ipt_script + \
          ' ' + ipt_obj_dir + \
          ' ' + ipt_path_to_preload_R + \
          ' ' + ipt_config_fl + \
          ' ' + output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)






if __name__ == "__main__":
    main()













