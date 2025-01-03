#!/usr/bin/env python

import argparse
import glob
import sys
import subprocess
import os

from input_other_required_scripts_dir.utils_motif_enrichment import s1_run_fimo
from input_other_required_scripts_dir.utils_motif_enrichment import s2_prepare_data


def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##step 01
    parser.add_argument("-genome_fl", dest='genome_fasta_file', help='Povide a genome fasta file.')

    parser.add_argument("-motif_fl", dest= 'motif_MEME_file',help = 'Provide a motif MEME file.')

    parser.add_argument("-acr_fl", dest = 'acr_file', help = 'Provide a ACR file used for predicting motifs within ACRs.')

    ##step 02
    parser.add_argument("-meta_fl", dest = 'meta_file', help = 'Provide a meta with annotated cell identity.')

    parser.add_argument("-tn5_fl", dest = 'tn5_file', help = 'Users will provide a BED file with the Tn5 integration sites.')

    parser.add_argument("-Gmfai_fl", dest='genome_fai_file', help='Provide a fai index file of reference genome.')

    ##step 02 and 03
    parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
                        help='Define the cluster column name of the meta file.')

    ##step 01 predict motif by fimo
    parser.add_argument("-open_predict_motif", dest='open_predict_motif_within_ACR', help='Perform motif prediction using fimo tool.'
                             'Default: yes')

    ##step 02 prepare the data for the enrichment
    parser.add_argument("-open_prepare_data", dest = 'open_prepare_data_enrichment', help = 'Prepare data before performing the enrichment test.'
                                                                                            'Default: yes')

    ##step 03 perform the enrichment test
    parser.add_argument("-open_enrichment", dest = 'open_enrichment_test', help = 'Perform the motif enrichment test.'
                                                                                  'Default: yes')

    ##optional parameters
    ##step 01
    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument("-qval", dest='q_value_cutoff_fimo', help = 'q value cutoff for the fimo')

    parser.add_argument("-ext_size", dest = 'acr_extended_size',help = '"Extended the flanking regions of ACRs to facilitate'
                                                                   ' the prediction of motifs within the ACRs and their surrounding sequences.'
                                            'Default: 5.')

    ##step 02
    parser.add_argument("-prefix", dest = 'sample_prefix', help = "Provide a prefix for the output."
                                                                  "Default: output")

    parser.add_argument("-downsampleCell", dest = 'down_sample_cell', help = 'Downsample cells to ensure a relative equal and unbiased comparison across cell types.'
                                                                             'Default: yes')

    parser.add_argument("-cell_num_cutoff", dest = 'cell_number_cutoff', help = 'Define a cutoff for the number of cells when downsampling.'
                                                                                'Default: 800')





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

    ##step 01
    if args.open_predict_motif_within_ACR is None:
        open_predict_motif_within_ACR_final = 'yes'
    else:
        if args.open_predict_motif_within_ACR == 'yes':
            open_predict_motif_within_ACR_final = 'yes'
        else:
            if args.open_predict_motif_within_ACR == 'no':
                open_predict_motif_within_ACR_final = 'no'
            else:
                print("Please use \'-open_predict_motif no\', to close this step")
                return

    ##step 02
    if args.open_prepare_data_enrichment is None:
        open_prepare_data_enrichment_final = 'yes'
    else:
        if args.open_prepare_data_enrichment == 'yes':
            open_prepare_data_enrichment_final = 'yes'
        else:
            if args.open_prepare_data_enrichment == 'no':
                open_prepare_data_enrichment_final = 'no'
            else:
                print("Please use \'-open_prepare_data no\', to close this step")
                return

    if args.down_sample_cell is None:
        down_sample_cell_final = 'yes'
    else:
        if args.down_sample_cell == 'yes':
            down_sample_cell_final = 'yes'
        else:
            if args.down_sample_cell == 'no':
                down_sample_cell_final = 'no'
            else:
                print("Please use \'-downsampleCell no\', to close this step")
                return

    ##step 03
    if args.open_enrichment_test is None:
        open_enrichment_test_final = 'yes'
    else:
        if args.open_enrichment_test == 'yes':
            open_enrichment_test_final = 'yes'
        else:
            if args.open_enrichment_test == 'no':
                open_enrichment_test_final = 'no'
            else:
                print("Please use \'-open_enrichment no\', to close this step")
                return


    if open_predict_motif_within_ACR_final == 'yes':

        if args.genome_fasta_file is None:
            print('Cannot find genome fasta file, please provide it')
            return
        else:
            try:
                file = open(args.genome_fasta_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the genome fasta file!')
                return

        if args.motif_MEME_file is None:
            print('Cannot find motif MEME file, please provide it')
            return
        else:
            try:
                file = open(args.motif_MEME_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the motif MEME file!')
                return

        if args.acr_file is None:
            print('Cannot find ACR file, please provide it')
            return
        else:
            try:
                file = open(args.acr_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the ACR file!')
                return

    if open_prepare_data_enrichment_final == 'yes':

        if args.meta_file is None:
            print('Cannot find cell meta file, please provide it')
            return
        else:
            try:
                file = open(args.meta_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the cell meta file!')
                return

        if args.tn5_file is None:
            print('Cannot find tn5 file, please provide it')
            return
        else:
            try:
                file = open(args.tn5_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the tn5 file!')
                return

        if args.genome_fai_file is None:
            print('Cannot find genome fai index file, please provide it')
            return
        else:
            try:
                file = open(args.genome_fai_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the fai index file!')
                return

    if open_enrichment_test_final == 'yes':

        if args.celltype_cluster_col_name is None:
            print('Please provide the name of column specifying the cell type identity in the meta file.')


    if args.required_script_dir is None:
        print('Cannot find required script dir, please provide the dir in \'-script_dir\' !')
        return

    input_required_scripts_dir = args.required_script_dir

    if args.core_number is not None:
        core_number_run = args.core_number
    else:
        core_number_run = '1'

    if args.q_value_cutoff_fimo is not None:
        q_value_cutoff_fimo_final = args.q_value_cutoff_fimo
    else:
        q_value_cutoff_fimo_final = '0.00005'

    if args.acr_extended_size is not None:
        acr_extended_size_final = args.acr_extended_size
    else:
        acr_extended_size_final = '5'

    if args.sample_prefix is not None:
        sample_prefix_final = args.sample_prefix
    else:
        sample_prefix_final = 'output'

    if args.cell_number_cutoff is not None:
        cell_number_cutoff_final = args.cell_number_cutoff
    else:
        cell_number_cutoff_final = '800'


    if open_predict_motif_within_ACR_final == 'yes':

        print ('Predict motif within ACRs')

        open_predict_motif_within_ACR_final_dir = output_dir + '/open_predict_motif_within_ACR_final_dir'
        if not os.path.exists(open_predict_motif_within_ACR_final_dir):
            os.makedirs(open_predict_motif_within_ACR_final_dir)

        input_acr_fl = args.acr_file
        input_genome_fasta_fl = args.genome_fasta_file
        input_motif_fl = args.motif_MEME_file

        s1_run_fimo.pipeline_analyze_acr_motif(input_acr_fl, input_genome_fasta_fl, input_motif_fl, open_predict_motif_within_ACR_final_dir,
                                               core_number_run, acr_extended_size_final)

        s1_run_fimo.wrap_all_results (open_predict_motif_within_ACR_final_dir,q_value_cutoff_fimo_final)

    if open_prepare_data_enrichment_final == 'yes':

        print ('Prepare data before enrichment test')

        open_prepare_data_enrichment_final_dir = output_dir + '/open_prepare_data_enrichment_final_dir'
        if not os.path.exists(open_prepare_data_enrichment_final_dir):
            os.makedirs(open_prepare_data_enrichment_final_dir)

        ###############
        input_peak_fl = args.acr_file
        input_tn5_bed_fl = args.tn5_file
        input_genome_fai_fl = args.genome_fai_file
        input_fastSparsetn5_pl = input_required_scripts_dir + '/utils_motif_enrichment/fastSparse.tn5.pl'
        prefix = sample_prefix_final

        s2_prepare_data.prepare_peak_tn5_sparse (input_peak_fl,input_tn5_bed_fl,input_genome_fai_fl,input_fastSparsetn5_pl,
                                                 prefix,
                                                 open_prepare_data_enrichment_final_dir)

        #############
        input_meta_fl = args.meta_file
        input_peak_cell_sparse_fl = open_prepare_data_enrichment_final_dir + '/opt_peak_sparse_sorted_dir/opt_peak_' + prefix + '_sorted.sparse'
        prefix_meta = 'addACRnum'

        s2_prepare_data.build_peak_cell_in_meta (input_meta_fl,input_peak_cell_sparse_fl,open_prepare_data_enrichment_final_dir,prefix_meta)

        if down_sample_cell_final == 'yes':

            cutoff_cellnm = cell_number_cutoff_final
            target_colnm = args.celltype_cluster_col_name

            s2_prepare_data.downsampling (cutoff_cellnm,open_prepare_data_enrichment_final_dir,prefix_meta,target_colnm)

            ##final output is opt_downsample_meta.txt


        ###############
        ipt_R_script = input_required_scripts_dir + '/utils_motif_enrichment/s2_prepare_regression_enrich.R'
        input_motif_acr_fl = output_dir + '/open_predict_motif_within_ACR_final_dir' + '/opt_final_ACR_motif.bed'
        input_motif_flt_cutoff = q_value_cutoff_fimo_final

        s2_prepare_data.prepare_motif_acr_data (ipt_R_script,input_motif_acr_fl,input_motif_flt_cutoff,input_peak_cell_sparse_fl,prefix,open_prepare_data_enrichment_final_dir)

    if open_enrichment_test_final == 'yes':

        print ('Perform motif enrichment test')

        open_enrichment_test_final_dir = output_dir + '/open_enrichment_test_final_dir'
        if not os.path.exists(open_enrichment_test_final_dir):
            os.makedirs(open_enrichment_test_final_dir)

        ipt_R_script = input_required_scripts_dir + '/utils_motif_enrichment/s3_motif_enrichment_test.R'

        if down_sample_cell_final == 'yes':
            ipt_final_meta_fl = output_dir + '/open_prepare_data_enrichment_final_dir' + '/opt_downsample_meta.txt'
        else:
            prefix_meta = 'addACRnum'
            ipt_final_meta_fl = output_dir + '/open_prepare_data_enrichment_final_dir' + '/opt_' + prefix_meta + '.txt'

        ipt_target_cluster_nm = args.celltype_cluster_col_name

        cmd = 'Rscript ' + ipt_R_script + \
              ' ' + output_dir + '/open_prepare_data_enrichment_final_dir' + \
              ' ' + ipt_final_meta_fl + \
              ' ' + open_enrichment_test_final_dir + \
              ' ' + ipt_target_cluster_nm + \
              ' ' + core_number_run
        print(cmd)
        subprocess.call(cmd, shell=True)





if __name__ == "__main__":
    main()


















