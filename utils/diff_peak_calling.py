#!/usr/bin/env python


##this script is to perform the peak calling
##it has two step
##step01 prepare the peak sparse file
##step02 make the DAR analysis
##updating 051925 step03 find the diff peak among diff groups
##updating 063025 we will use the object to


import argparse
import glob
import sys
import os
import subprocess
import re


from input_other_required_scripts_dir.utils_diff_peak_calling_python import s1_make_peak_sparse as s1_prepare


def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    ##step 01
    parser.add_argument("-s1_open_prepare", dest='s1_open_prepare_peak_tn5',
                        help='Run the step 01 to prepare peak tn5 file.'
                             'Default: yes')

    ##updating 063025
    parser.add_argument("-soc_obj", dest='soc_object_fl',
                        help='Provide an object obtained from the peak calling step.')

    parser.add_argument("-peak_fl", dest = 'peak_file', help = 'Provide the final peak file.')

    parser.add_argument("-tn5_fl", dest='tn5_bed_file', help="Provide the tn5 bed file.")


    ##step 02
    parser.add_argument("-s2_open_diff_ct_peak", dest='s2_open_diff_peak_call',help='Run the step 02 to call the DA peak per cell type.'
                                                                                 'Default: yes')


    parser.add_argument("-meta_fl", dest='meta_file', help='Provide the meta file recording cell identity information.')


    ##updating 051925
    ##step 03
    parser.add_argument("-s3_open_diff_gp_peak", dest= 's3_open_diff_group_peak', help = 'Run the step 03 to call the DA peak for the specific group.'
                                                                                         'Default: no')


    ##Optional
    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
                        help='Define the cluster column name of the meta file.')

    parser.add_argument("-stat_test", dest='stat_test_method', help='Users need to provide the method used for the stat test. pnorm or perm.'
                                                                    'Default: perm.')

    parser.add_argument("-threshold", dest='threshold_val', help = 'The threshold value for statistical testing.'
                                                                   'Default: 0.001')

    parser.add_argument("-null_permutations", dest='null_permutations_val', help = 'The number of null permutations to generate.'
                                                                   'Default: 1000')

    parser.add_argument("-entropy_bootstraps", dest = 'entropy_bootstraps_val', help = 'The number of bootstraps for the entropy metric to generate.'
                                                                                       'Default: 1000')

    parser.add_argument("-prefix", dest = 'prefix_str', help = 'Prefix for output.'
                                                               'Default:opt')


    ##updating 051925
    parser.add_argument("-gp_colnm", dest = 'group_col_name', help = 'Users need to provide a group column name will be used for the diff peaks calling.')


    #parser.add_argument("-upsample_method", dest= 'upsampling_method',help = 'Users need to define the method used for the upsampling. Sample read or sample cell'
    #                                               'Default: read')

    #parser.add_argument("-FDR", dest = 'FDR_cutoff', help = 'Provide a FDR cutoff to filter the differentially accessible peaks'
    #                                                        'Default: 0.1')

    #parser.add_argument("-log2fc", dest= 'log2fc_cutoff', help = 'Provde a log2fc cutoff to filter the differentially accessible peaks.'
    #                                                             'Default: 1')

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
    else:
        input_required_scripts_dir = args.required_script_dir


    if args.s1_open_prepare_peak_tn5 is None:

        s1_open_prepare_peak_tn5_final = 'yes'

    else:

        if args.s1_open_prepare_peak_tn5 == 'yes':

            s1_open_prepare_peak_tn5_final = 'yes'

            ##updating 063025
            ##if users provide the soc from the last step
            if args.soc_object_fl is not None:
                print('Use object from the last step as the input file')
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the object file!')
                    return
                use_soc_object_as_ipt = 'yes'
            else:
                use_soc_object_as_ipt = 'no'

            if use_soc_object_as_ipt == 'no':

                print('Do not use object from the last step as the input file and please provide the peak file as input')

                if args.peak_file is None:
                    print('Cannot find the peak file, please provide it')
                    return
                else:
                    try:
                        file = open(args.peak_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the peak file!')
                        return

            if args.tn5_bed_file is None:
                print('Cannot find the tn5 file, please provide it')
                return
            else:
                try:
                    file = open(args.tn5_bed_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the tn5 file!')
                    return

        else:
            s1_open_prepare_peak_tn5_final = 'no'
            print(
                'Users choose to close the peak preparation, please use \'-s1_open_prepare yes\' to open this step')


    if args.s2_open_diff_peak_call is None:

        s2_open_diff_peak_call_final = 'yes'

    else:

        if args.s2_open_diff_peak_call == 'yes':

            s2_open_diff_peak_call_final = 'yes'

            ##updating 063025
            if args.soc_object_fl is not None:
                print('Use object from the last step as the input file')
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the object file!')
                    return
                use_soc_object_as_ipt = 'yes'
            else:
                use_soc_object_as_ipt = 'no'

            if use_soc_object_as_ipt == 'no':

                print(
                    'Do not use object from the last step as the input file and please provide the meta file as input')

                if args.meta_file is None:
                    print('Cannot find the meta file, please provide it')
                    return
                else:
                    try:
                        file = open(args.meta_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the meta file!')
                        return


        else:

            s2_open_diff_peak_call_final = 'no'
            print('Users choose to close the open diff peak calling, please use \'-s2_open_diff_peak yes\' to open this step')


    ##updating 051925
    if args.s3_open_diff_group_peak is None:

        s3_open_diff_group_peak_final = 'no'

    else:

        if args.s3_open_diff_group_peak == 'yes':
            s3_open_diff_group_peak_final = 'yes'

            if args.group_col_name is None:
                print('Please provide a group column name within the meta file')
                return


        else:
            print('Users choose to close the open diff peak calling per group, please use \'-s3_open_diff_gp_peak yes\' to open this step')
            return



    store_final_parameter_line_list = []
    ##step02 parameters
    if args.core_number is not None:
        core_number_final = args.core_number
    else:
        core_number_final = '1'

    store_final_parameter_line_list.append('thread_num <- ' + core_number_final)


    if args.celltype_cluster_col_name is not None:
        celltype_cluster_col_name_final = args.celltype_cluster_col_name
        store_final_parameter_line_list.append('target_cluster <- ' + '\'' + celltype_cluster_col_name_final + '\'')
    else:

        ##updating 063025
        if args.soc_object_fl is not None:
            celltype_cluster_col_name_final = 'cell_identity'
            store_final_parameter_line_list.append('target_cluster <- ' + '\'' + celltype_cluster_col_name_final + '\'')
        else:
            if s2_open_diff_peak_call_final == 'yes':
                print ('Please use \'-ct_colnm\' to specify the column name showing cell identity once users do not use the object as the input file.')
                return


    if args.stat_test_method is not None:
        stat_test_method_final = args.stat_test_method
    else:
        stat_test_method_final = 'perm'

    store_final_parameter_line_list.append('stat_test <- ' + '\'' + stat_test_method_final + '\'')


    if args.threshold_val is not None:
        threshold_val_final = args.threshold_val
    else:
        threshold_val_final = '0.0001'

    store_final_parameter_line_list.append('threshold <- ' + threshold_val_final)


    if args.null_permutations_val is not None:
        null_permutations_val_final = args.null_permutations_val
    else:
        null_permutations_val_final = '1000'

    store_final_parameter_line_list.append('null_permutations <- ' + null_permutations_val_final)


    if args.entropy_bootstraps_val is not None:
        entropy_bootstraps_val_final = args.entropy_bootstraps_val
    else:
        entropy_bootstraps_val_final = '1000'

    store_final_parameter_line_list.append('entropy_bootstraps <- ' + entropy_bootstraps_val_final)


    if args.prefix_str is not None:
        prefix_str_final = args.prefix_str
    else:
        prefix_str_final = 'opt'

    store_final_parameter_line_list.append('prefix <- ' + '\'' + prefix_str_final + '\'')


    ##updating 051925
    if args.group_col_name is not None:
        group_col_name_final = args.group_col_name
    else:
        group_col_name_final = 'na'

    store_final_parameter_line_list.append('target_treat_colnm <- ' + '\'' + group_col_name_final + '\'')


    with open(output_dir + '/temp_defined_parameters.config', 'w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')


    if s1_open_prepare_peak_tn5_final == 'yes':

        print ('Users will prepare the peak tn5 sparse file before calling cell type specific peaks.')

        s1_open_prepare_peak_tn5_final_dir = output_dir + '/s1_open_prepare_peak_tn5_final'
        if not os.path.exists(s1_open_prepare_peak_tn5_final_dir):
            os.makedirs(s1_open_prepare_peak_tn5_final_dir)

        ##updating 063025
        if args.soc_object_fl is not None:

            ipt_object_file = args.soc_object_fl
            cmd = 'Rscript ' + input_required_scripts_dir + '/utils_diff_peak_calling_python/read_object.R' + \
                  ' ' + ipt_object_file + \
                  ' ' + s1_open_prepare_peak_tn5_final_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

            input_peak_fl = s1_open_prepare_peak_tn5_final_dir + '/temp_peak.txt'

        else:

            input_peak_fl = args.peak_file

        input_tn5_bed_fl = args.tn5_bed_file
        input_fastSparsetn5_pl = input_required_scripts_dir + '/utils_diff_peak_calling_python/fastSparse.tn5.py'

        s1_prepare.make_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
                                s1_open_prepare_peak_tn5_final_dir)

        s1_prepare.sort_sparse(s1_open_prepare_peak_tn5_final_dir)

        input_peak_tn5_fl = s1_open_prepare_peak_tn5_final_dir + '/opt_peak_sparse_sorted_dir/opt_peak_sorted.sparse'

        s1_prepare.prepare_peak_acc (input_peak_tn5_fl,input_peak_fl,s1_open_prepare_peak_tn5_final_dir)

        ##final output for the next step is storing
        ##opt_prepare_peak_acc_dir/

    if s2_open_diff_peak_call_final == 'yes':

        print ('Users will call the differential accessible peaks across cell types.')

        s2_open_diff_peak_call_final_dir = output_dir + '/s2_open_diff_peak_call_final_dir'
        if not os.path.exists(s2_open_diff_peak_call_final_dir):
            os.makedirs(s2_open_diff_peak_call_final_dir)

        R_script = input_required_scripts_dir + '/utils_diff_peak_calling_python/s2_call_ctACR.Bootstrap.R'
        R_script_saveobj = input_required_scripts_dir + '/utils_diff_peak_calling_python/s2_call_ctACR.Bootstrap_saveobj.R'

        ipt_peak_sparse_fl = output_dir + '/s1_open_prepare_peak_tn5_final/opt_prepare_peak_acc_dir/opt_accessibility.txt'

        ##updating 063025
        if args.soc_object_fl is not None:
            ipt_meta_fl = output_dir + '/s1_open_prepare_peak_tn5_final/temp_unmodi_update_meta.txt'
        else:
            ipt_meta_fl = args.meta_file

        ipt_peak_fl = output_dir + '/s1_open_prepare_peak_tn5_final/opt_prepare_peak_acc_dir/opt_peak.bed'

        if os.path.isfile(ipt_peak_sparse_fl) == True:

            ##updating 063025
            if args.soc_object_fl is not None:

                ipt_soc_obj_fl = args.soc_object_fl

                if re.match('.+/(.+)\.atac\.soc\.rds', ipt_soc_obj_fl):
                    mt = re.match('.+/(.+)\.atac\.soc\.rds', ipt_soc_obj_fl)
                    input_prefix = mt.group(1)
                else:
                    if re.match('(.+)\.atac\.soc\.rds', ipt_soc_obj_fl):
                        mt = re.match('(.+)\.atac\.soc\.rds', ipt_soc_obj_fl)
                        input_prefix = mt.group(1)
                    else:
                        print('Please use *.atac.soc.rds file without changing the file name')
                        return

                cmd = 'Rscript ' + R_script_saveobj + \
                      ' ' + ipt_peak_sparse_fl + \
                      ' ' + ipt_meta_fl + \
                      ' ' + ipt_peak_fl + \
                      ' ' + ipt_soc_obj_fl + \
                      ' ' + input_prefix + \
                      ' ' + output_dir + '/temp_defined_parameters.config' + \
                      ' ' + s2_open_diff_peak_call_final_dir
                print(cmd)
                subprocess.call(cmd,shell=True)


            else:
                ##check the peak sparse fl is existing or not
                cmd = 'Rscript ' + R_script + \
                      ' ' + ipt_peak_sparse_fl + \
                      ' ' + ipt_meta_fl + \
                      ' ' + ipt_peak_fl + \
                      ' ' + output_dir + '/temp_defined_parameters.config' + \
                      ' ' + s2_open_diff_peak_call_final_dir
                print(cmd)
                subprocess.call(cmd,shell=True)

        else:
            print('Please check the step01 and make sure opt_accessibility.txt and opt_peak.bed file existing in the s1_open_prepare_peak_tn5_final_dir/opt_prepare_peak_acc_dir')



    if s3_open_diff_group_peak_final == 'yes':

        print('Users will call the differential accessible peaks across defined group.')

        s3_open_diff_group_peak_final_dir = output_dir + '/s3_open_diff_group_peak_final_dir'
        if not os.path.exists(s3_open_diff_group_peak_final_dir):
            os.makedirs(s3_open_diff_group_peak_final_dir)

        R_script = input_required_scripts_dir + '/utils_diff_peak_calling_python/s3_call_treat_ctACR.Bootstrap.R'

        ipt_peak_sparse_fl = output_dir + '/s1_open_prepare_peak_tn5_final/opt_prepare_peak_acc_dir/opt_accessibility.txt'

        ##updating 063025

        if args.soc_object_fl is not None:
            ipt_meta_fl = output_dir + '/s1_open_prepare_peak_tn5_final/temp_unmodi_update_meta.txt'
        else:
            ipt_meta_fl = args.meta_file


        ipt_peak_fl = output_dir + '/s1_open_prepare_peak_tn5_final/opt_prepare_peak_acc_dir/opt_peak.bed'

        if os.path.isfile(ipt_peak_sparse_fl) == True:

            cmd = 'Rscript ' + R_script + \
                  ' ' + ipt_peak_sparse_fl + \
                  ' ' + ipt_meta_fl + \
                  ' ' + ipt_peak_fl + \
                  ' ' + output_dir + '/temp_defined_parameters.config' + \
                  ' ' + s3_open_diff_group_peak_final_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

        else:
            print('Please check the step01 and make sure opt_accessibility.txt and opt_peak.bed file existing in the s1_open_prepare_peak_tn5_final_dir/opt_prepare_peak_acc_dir')






if __name__ == "__main__":
    main()







