#!/usr/bin/env python


##this is python version of calling peaks
##this script is to directly run MACS2 based on previous developed bed files
##and save the dtemp files to run MACS2 for different threshold

##updating 062925 we will change the meta and ct colnm into the obj
##updating 062425 we will filter peaks excess the range of Gmfai
##updating 110524 we will introduce a new function to prepare the big wig file

import argparse
import glob
import sys
import os
import subprocess
import re


from input_other_required_scripts_dir.utils_cell_type_peak_calling_dir_python import step01_prepare_tn5_file_per_celltype as s1_tn5
from input_other_required_scripts_dir.utils_cell_type_peak_calling_dir_python import step02_call_peak_per_celltype as s2_pcall
from input_other_required_scripts_dir.utils_cell_type_peak_calling_dir_python import step03_filter_combine_all_peak as s3_fpeak
from input_other_required_scripts_dir.utils_cell_type_peak_calling_dir_python import step04_build_bigwig_file as s4_bigwig

def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument("-tn5_fl", dest='tn5_bed_file', help="Provide the tn5 bed file.")

    parser.add_argument("-meta_fl", dest='meta_file', help='Provide the meta file recording cell identity information.'
                                                           'Make sure last column is the cell identity information')

    ##updating 062925
    parser.add_argument("-soc_obj", dest='soc_object_fl',
                        help='Provide an object obtained from the cell type annotation step.')


    parser.add_argument("-Ggff_fl", dest='gene_gff_file', help='Provide the gene gff file used to extract exon region when permutating peaks.')

    parser.add_argument("-Gmfai_fl", dest = 'genome_fai_file', help = 'Provide a fai index file of reference genome.')

    parser.add_argument("-script_dir", dest = 'required_script_dir', help = 'Users must provide the required script dir provided by GitHub.')

    #parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
    #                    help='Define the cluster column name of the meta file.')

    ##optional parameters
    parser.add_argument("-s1_open_tn5", dest = 's1_open_prepare_tn5_file', help = 'Run the step 01 to prepare the tn5 file per cell type.'
                                                                                 'Default: yes')

    parser.add_argument("-s2_open_peak", dest = 's2_open_call_peak', help = 'Run the step 02 to call peaks per cell type.'
                                                                           'Default: yes')

    parser.add_argument("-s3_open_filter", dest = 's3_open_filter_peak', help = 'Run the step 03 to filter peaks from all the peak calling.'
                                                                               'Default: yes')

    parser.add_argument("-s4_open_bigwig", dest = 's4_open_bigwig_build', help = 'Run the step 04 to obtain the bigwig file per cell type.'
                                                                                 'Default: yes')

    ##updating 062925
    parser.add_argument("-open_build_object", dest = 'open_build_object_file', help = 'After deciding the final peak file, please use this step to build the final object.'
                                                                                      'Default: no')

    parser.add_argument("-final_peak_fl" , dest = 'final_peak_file', help = 'Specify the final peak file after deciding which peak file is suitable for the peak calling.')



    parser.add_argument("-core", dest = 'core_number', help = 'Specify how many cores we will use.'
                                                              'Default: 1')

    ##updating 062925
    parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
                        help='Define the cluster column name of the meta file. Default: cell_identity.')


    ##step02 parameters
    parser.add_argument("-umapp_fl", dest = 'unmappable_file', help = 'Provide the unmappable file that will help build the imputated region.'
                                                                      'Default: None')

    parser.add_argument("-g_size", dest = 'genome_size', help = 'Provide an effective genome size. Users must provide it.')

    parser.add_argument("-p_or_q", dest='pval_or_qval', help= 'Provide to use pvalue or qvalue cutoff using MACS2.'
                                                              'Default: pval')

    parser.add_argument("-cutoff", dest='pqval_cutoff', help='Provide an exact cutoff for either pvalue or qvalue.'
                                                             'Default: 0.05')

    parser.add_argument("-SLOCAL", dest='SLOCAL_val', help='Provide a parameter setting of SLOCAL within MACS2 to adjust peak calling sensistivity.'
                                                           'Default: 75')

    parser.add_argument("-LLOCAL", dest='LLOCAL_val', help='Provide a parameter setting of LLOCAL within MACS2 to adjust peak calling sensistivity.'
                                                           'Default: 500')

    parser.add_argument("-max_gap", dest='max_gap_val', help='Provide a parameter setting of max_gap within MACS2 to adjust peak calling sensistivity.'
                                                           'Default: 10')

    parser.add_argument("-FDR", dest = 'FDR_string', help = 'Provide a FDR string to return all the peak files with all defined FDR values.'
                                                            'Default: 0.05,0.01,0,005,0.001')



    ##parse of parameters
    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()



    if args.required_script_dir is None:
        print('Cannot find required script dir, please provide the dir in \'-script_dir\' !')
        return
    else:
        input_required_scripts_dir = args.required_script_dir


    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir


    ##udpating 062925
    open_use_soc_obj = 'no'
    open_use_meta_fl = 'no'
    if args.soc_object_fl is None:

        print('Cannot find the soc object file, please provide it. Or provide a meta file')

        if args.meta_file is None:
            print('Cannot find meta file, please provide it')
            return
        else:
            open_use_meta_fl = 'yes'
            try:
                file = open(args.meta_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the meta file!')
                return
    else:
        open_use_soc_obj = 'yes'
        try:
            file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the soc object file!')
            return


    if args.s1_open_prepare_tn5_file is None:
        s1_open_prepare_tn5_file_final = 'yes'
    else:
        if args.s1_open_prepare_tn5_file == 'yes':
            s1_open_prepare_tn5_file_final = 'yes'
        else:
            s1_open_prepare_tn5_file_final = 'no'
            print('Users close the step01, please use \'-s1_open_tn5\' yes to open this step')

    ##check if we provide the target cluster name
    #if s1_open_prepare_tn5_file_final == 'yes':

    #    if args.celltype_cluster_col_name is None:
    #        print('Please provide the name of column specifying the cell type identity in the meta file.')
    #        return

    if args.s2_open_call_peak is None:
        s2_open_call_peak_final = 'yes'
    else:
        if args.s2_open_call_peak == 'yes':
            s2_open_call_peak_final = 'yes'
        else:
            s2_open_call_peak_final = 'no'
            print('Users close the step02, please use \'-s2_open_peak\' yes to open this step')

    if args.s3_open_filter_peak is None:
        s3_open_filter_peak_final = 'yes'
    else:
        if args.s3_open_filter_peak == 'yes':
            s3_open_filter_peak_final = 'yes'
        else:
            s3_open_filter_peak_final = 'no'
            print('Users close the step03, please use \'-s3_open_filter\' yes to open this step')

    if args.s4_open_bigwig_build is None:
        s4_open_bigwig_build_final = 'yes'
    else:
        if args.s4_open_bigwig_build == 'yes':
            s4_open_bigwig_build_final = 'yes'
        else:
            s4_open_bigwig_build_final = 'no'
            print('Users close the step04, please use \'-s4_open_bigwig\' yes to open this step')


    if args.open_build_object_file is not None:

        if args.open_build_object_file != 'yes':
            print('Please specify the -open_build_object yes')
            return
        else:
            open_build_object_file_final = 'yes'

            if args.final_peak_file is None:
                print('Cannot find final peak file, please provide it')
                return
            else:
                try:
                    file = open(args.final_peak_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the final peak file!')
                    return

            if args.soc_object_fl is None:
                print('Cannot find soc object file, please provide it')
                return
            else:
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the soc object peak file!')
                    return

    else:
        open_build_object_file_final = 'no'


    if args.unmappable_file is not None:
        ##check the file
        input_unmappable_fl = args.unmappable_file
        try:
            file = open(args.unmappable_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the domain_file!')
            return
    else:
        input_unmappable_fl = 'None'


    if args.core_number is not None:
        core_number_run = args.core_number
    else:
        core_number_run = '1'




    if args.pval_or_qval is None:
        final_pval_or_qval = 'pval'
    else:
        final_pval_or_qval = args.pval_or_qval

    if args.pqval_cutoff is None:
        final_pqval_cutoff = '0.05'
    else:
        final_pqval_cutoff = args.pqval_cutoff

    if args.SLOCAL_val is None:
        final_SLOCAL_val = '75'
    else:
        final_SLOCAL_val = args.SLOCAL_val

    if args.LLOCAL_val is None:
        final_LLOCAL_val = '500'
    else:
        final_LLOCAL_val = args.LLOCAL_val

    if args.max_gap_val is None:
        final_max_gap_val = '10'
    else:
        final_max_gap_val = args.max_gap_val

    if args.FDR_string is None:
        final_FDR_string = '0.05,0.01,0.005,0.001'
    else:
        final_FDR_string = args.FDR_string

    ##updating 062925
    if args.celltype_cluster_col_name is None:
        target_colnm = 'cell_identity'
    else:
        target_colnm = args.celltype_cluster_col_name


    utils_cell_type_peak_calling_dir = input_required_scripts_dir + '/utils_cell_type_peak_calling_dir_python'


    ##updating 062925
    if open_build_object_file_final == 'no':

        ##check the file
        if args.genome_fai_file is None:
            print('Cannot find genome fai index file, please provide it')
            return
        else:
            try:
                file = open(args.genome_fai_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the genome fai index file!')
                return


        if s1_open_prepare_tn5_file_final == 'yes':

            ##check the file
            if args.tn5_bed_file is None:
                print('Cannot find tn5 file, please provide it')
                return
            else:
                try:
                    file = open(args.tn5_bed_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the tn5 file!')
                    return


            print('Users choose to run the step 01 to prepare the tn5 bed file per cell type')

            s1_open_prepare_tn5_file_final_dir = output_dir + '/s1_open_prepare_tn5_file_final_dir'
            if not os.path.exists(s1_open_prepare_tn5_file_final_dir):
                os.makedirs(s1_open_prepare_tn5_file_final_dir)

            input_tn5_fl = args.tn5_bed_file
            #input_meta_fl = args.meta_file
            input_core_num = core_number_run
            #target_colnm = args.celltype_cluster_col_name

            ##updating 062925
            if open_use_soc_obj == 'yes':
                input_soc_obj_fl = args.soc_object_fl
                cmd = 'Rscript ' + utils_cell_type_peak_calling_dir + '/read_object.R' + \
                      ' ' + input_soc_obj_fl + \
                      ' ' + s1_open_prepare_tn5_file_final_dir
                print(cmd)
                subprocess.call(cmd,shell=True)
                ipt_meta_fl = s1_open_prepare_tn5_file_final_dir + '/temp_unmodi_update_meta.txt'

            else:
                if open_use_meta_fl == 'yes':
                    ipt_meta_fl = args.meta_file
                else:
                    print('Please specify the -soc_obj yes or -meta_fl yes to make a meta file')
                    return

            ##updating 031525
            ##modify the meta file to allow the cluster cell type name without ' '
            store_final_line_list = []
            count = 0
            with open(ipt_meta_fl, 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    count += 1
                    if count == 1:
                        store_final_line_list.append(eachline)
                    else:
                        eachline = eachline.replace(' ','_')
                        store_final_line_list.append(eachline)

            with open (s1_open_prepare_tn5_file_final_dir + '/temp_update_meta.txt','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            input_meta_fl = s1_open_prepare_tn5_file_final_dir + '/temp_update_meta.txt'

            s1_tn5.prepare_bed_files_parallele(input_tn5_fl, input_meta_fl,target_colnm,
                                        input_core_num,
                                        s1_open_prepare_tn5_file_final_dir)

        if s2_open_call_peak_final == 'yes':

            ##check the file

            if args.gene_gff_file is None:
                print('Cannot find gene gff file, please provide it')
                return
            else:
                try:
                    file = open(args.gene_gff_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the gene gff file!')
                    return

            if args.genome_size is None:
                print('Cannot find effect genome size argument, please provide it')
                return
            else:
                g_size = args.genome_size

            print('Users choose to run the step 02 to run the peak calling per cell type')

            s2_open_call_peak_final_dir = output_dir + '/s2_open_call_peak_final_dir'
            if not os.path.exists(s2_open_call_peak_final_dir):
                os.makedirs(s2_open_call_peak_final_dir)

            input_pool_bed_dir = output_dir + '/s1_open_prepare_tn5_file_final_dir/store_pool_bed_dir'

            all_temp_output_dir_list = glob.glob(input_pool_bed_dir + '/*')
            ipt_target_bed_fl_list = []
            for eachdir in all_temp_output_dir_list:
                all_bed_fl_list = glob.glob(eachdir + '/*')
                for eachbedfl in all_bed_fl_list:
                    ipt_target_bed_fl_list.append(eachbedfl)

            input_gff_fl = args.gene_gff_file
            input_ref_fai_fl = args.genome_fai_file
            input_core_num = core_number_run
            shiftsize = '-75'
            extsize = '150'
            open_callpeak = 'yes'
            open_permute = 'yes'

            s2_pcall.run_pipeline(s2_open_call_peak_final_dir, ipt_target_bed_fl_list, input_gff_fl, input_ref_fai_fl,
                                 input_unmappable_fl,
                                 input_core_num,
                                 g_size,
                                 final_pval_or_qval,
                                 final_pqval_cutoff,
                                 final_SLOCAL_val, final_LLOCAL_val, final_max_gap_val, extsize, shiftsize,
                                 open_callpeak, open_permute)

        if s3_open_filter_peak_final == 'yes':

            print('Users choose to filter and merge peaks based on the previous peak calling')

            s3_open_filter_peak_final_dir = output_dir + '/s3_open_filter_peak_final_dir'
            if not os.path.exists(s3_open_filter_peak_final_dir):
                os.makedirs(s3_open_filter_peak_final_dir)

            s2_targettn5_colNum = '2'
            #input_meta_fl = args.meta_file
            #target_colnm = args.celltype_cluster_col_name

            input_meta_fl = output_dir + '/s1_open_prepare_tn5_file_final_dir' + '/temp_update_meta.txt'

            finalfdr_list = final_FDR_string.split(',')

            s3_fpeak.FDR_filteration(output_dir + '/s2_open_call_peak_final_dir', utils_cell_type_peak_calling_dir,
                            input_meta_fl, finalfdr_list, s2_targettn5_colNum,s3_open_filter_peak_final_dir,final_pval_or_qval,final_pqval_cutoff,target_colnm)

            ##updating 062425
            input_ref_fai_fl = args.genome_fai_file
            all_fl_list = glob.glob(s3_open_filter_peak_final_dir + '/*Peaks.bed')
            for eachfl in all_fl_list:
                s3_fpeak.filter_peak_by_chr_size(eachfl, input_ref_fai_fl, s3_open_filter_peak_final_dir)



        if s4_open_bigwig_build_final == 'yes':

            print('Users choose to build up the bigwig file')

            s4_open_bigwig_build_final = output_dir + '/s4_open_bigwig_build_final_dir'
            if not os.path.exists(s4_open_bigwig_build_final):
                os.makedirs(s4_open_bigwig_build_final)

            s4_norm_type = 'tn5count'
            input_ref_fai_fl = args.genome_fai_file
            input_core_num = core_number_run
            s4_open_sort_bdg = 'no'

            s4_bigwig.build_bigwig_fl_parallele(output_dir,s4_open_bigwig_build_final, utils_cell_type_peak_calling_dir, input_ref_fai_fl, input_core_num,
                                       s4_open_sort_bdg, s4_norm_type)


    else:

        print('Users choose to run the step to save the object')

        ipt_final_peak_fl = args.final_peak_file
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

        cmd = 'Rscript ' + utils_cell_type_peak_calling_dir + '/save_object.R' + \
              ' ' + ipt_soc_obj_fl + \
              ' ' + ipt_final_peak_fl + \
              ' ' + output_dir + \
              ' ' + input_prefix
        print(cmd)
        subprocess.call(cmd,shell=True)




if __name__ == "__main__":
    main()






















