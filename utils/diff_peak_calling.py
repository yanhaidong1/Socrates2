#!/usr/bin/env python


##this script is to perform the peak calling
##it has two step
##step01 prepare the peak sparse file
##step02 make the DAR analysis

import argparse
import glob
import sys
import os
import subprocess


from input_other_required_scripts_dir.utils_diff_peak_calling import s1_make_peak_sparse as s1_prepare


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

    parser.add_argument("-peak_fl", dest = 'peak_file', help = 'Provide the final peak file.')

    parser.add_argument("-tn5_fl", dest='tn5_bed_file', help="Provide the tn5 bed file.")


    ##step 02
    parser.add_argument("-s2_open_diff_peak", dest='s2_open_diff_peak_call',help='Run the step 02 to call the DA peak.'
                                                                                 'Default: yes')


    parser.add_argument("-meta_fl", dest='meta_file', help='Provide the meta file recording cell identity information.')


    ##Optional
    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
                        help='Define the cluster column name of the meta file.')

    parser.add_argument("-upsample_method", dest= 'upsampling_method',help = 'Users need to define the method used for the upsampling. Sample read or sample cell'
                                                   'Default: read')

    parser.add_argument("-FDR", dest = 'FDR_cutoff', help = 'Provide a FDR cutoff to filter the differentially accessible peaks'
                                                            'Default: 0.1')

    parser.add_argument("-log2fc", dest= 'log2fc_cutoff', help = 'Provde a log2fc cutoff to filter the differentially accessible peaks.'
                                                                 'Default: 1')

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

        s1_open_prepare_peak_tn5_final = 'no'

    else:

        if args.s1_open_prepare_peak_tn5 == 'yes':

            s1_open_prepare_peak_tn5_final = 'yes'

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

        s2_open_diff_peak_call_final = 'no'

    else:

        if args.s2_open_diff_peak_call == 'yes':

            s2_open_diff_peak_call_final = 'yes'

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
        if s2_open_diff_peak_call_final == 'yes':
            print ('Please use \'-ct_colnm\' to specify the column name showing cell identity.')
            return


    if args.upsampling_method is not None:
        upsampling_method_final = args.upsampling_method
        if upsampling_method_final != 'read' and upsampling_method_final != 'cell':
            print ('Please make sure \'-upsample_method\' should be \'read\' or \'cell\'')
            return
    else:
        upsampling_method_final = 'read'

    store_final_parameter_line_list.append('select_method <- ' + '\'' + upsampling_method_final + '\'')


    if args.FDR_cutoff is not None:
        FDR_cutoff_final = args.FDR_cutoff
    else:
        FDR_cutoff_final = '0.1'

    store_final_parameter_line_list.append('FDR_cutoff <- ' + FDR_cutoff_final)

    if args.log2fc_cutoff is not None:
        log2fc_cutoff_final = args.log2fc_cutoff
    else:
        log2fc_cutoff_final = '1'

    store_final_parameter_line_list.append('log2fc_cutoff <- ' + log2fc_cutoff_final)

    with open(output_dir + '/temp_defined_parameters.config', 'w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')


    if s1_open_prepare_peak_tn5_final == 'yes':

        print ('Users will prepare the peak tn5 sparse file before calling cell type specific peaks.')

        s1_open_prepare_peak_tn5_final_dir = output_dir + '/s1_open_prepare_peak_tn5_final'
        if not os.path.exists(s1_open_prepare_peak_tn5_final_dir):
            os.makedirs(s1_open_prepare_peak_tn5_final_dir)


        input_peak_fl = args.peak_file
        input_tn5_bed_fl = args.tn5_bed_file
        input_fastSparsetn5_pl = input_required_scripts_dir + '/utils_diff_peak_calling/fastSparse.tn5.pl'

        s1_prepare.make_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
                                s1_open_prepare_peak_tn5_final_dir)

        s1_prepare.sort_sparse(s1_open_prepare_peak_tn5_final_dir)

    if s2_open_diff_peak_call_final == 'yes':

        print ('Users will call the differential accessible peaks across cell types.')

        s2_open_diff_peak_call_final_dir = output_dir + '/s2_open_diff_peak_call_final_dir'
        if not os.path.exists(s2_open_diff_peak_call_final_dir):
            os.makedirs(s2_open_diff_peak_call_final_dir)

        R_script = input_required_scripts_dir + '/utils_diff_peak_calling/s2_pseudobulk_DAR_analysis.R'
        ipt_peak_sparse_fl = output_dir + '/s1_open_prepare_peak_tn5_final/opt_peak_sparse_sorted_dir/opt_peak_sorted.sparse'
        ipt_meta_fl = args.meta_file

        if os.path.isfile(ipt_peak_sparse_fl) == True:

            ##check the peak sparse fl is existing or not
            cmd = 'Rscript ' + R_script + \
                  ' ' + ipt_peak_sparse_fl + \
                  ' ' + ipt_meta_fl + \
                  ' ' + s2_open_diff_peak_call_final_dir + \
                  ' ' + output_dir + '/temp_defined_parameters.config'
            print(cmd)
            subprocess.call(cmd,shell=True)

        else:
            print('Please check the step01 and make sure opt_peak_sorted.sparse file exists in the s1_open_prepare_peak_tn5_final_dir/opt_peak_sparse_sorted_dir')





if __name__ == "__main__":
    main()







