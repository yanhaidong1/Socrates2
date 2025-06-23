#!/usr/bin/env python

import argparse
import glob
import sys
import subprocess
import os
from Bio import SeqIO

from input_other_required_scripts_dir.utils_motif_analysis import s1_run_fimo
from input_other_required_scripts_dir.utils_motif_analysis import s2_prepare_data


##updating 040825 we will add a new function to perform motif deviation analysis and perform the motif UMAP plotting


def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##updating 040925
    ##add a step to build the chromVAR
    parser.add_argument("-open_motif_enrich", dest='open_motif_enrichment_analysis', help = 'Open the motif enrichment analysis to obtain the motifs enriched in the cell-type-specific ACRs.')

    parser.add_argument("-open_motif_DevScore", dest = 'open_motif_deviation_score', help = 'Open the calculation of motif deviation score.')

    parser.add_argument("-open_build_BSgenome", dest='open_build_BSgenome_lib',
                        help='This step helps users build a BSgenome library.')




    ##step 01
    parser.add_argument("-genome_fl", dest='genome_fasta_file', help='Povide a genome fasta file.')

    parser.add_argument("-motif_fl", dest= 'motif_MEME_file',help = 'Provide a motif MEME file.')

    parser.add_argument("-acr_fl", dest = 'acr_file', help = 'Provide a ACR file used for predicting motifs within ACRs.')

    ##step 02
    parser.add_argument("-meta_fl", dest = 'meta_file', help = 'Provide a meta with annotated cell identity.')

    parser.add_argument("-tn5_fl", dest = 'tn5_file', help = 'Users will provide a BED file with the Tn5 integration sites.')

    #parser.add_argument("-Gmfai_fl", dest='genome_fai_file', help='Provide a fai index file of reference genome.')

    ##step 02 and 03
    parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
                        help='Define the cluster column name of the meta file.')

    ##step 01 predict motif by fimo
    parser.add_argument("-s1_open_predict_motif", dest='open_predict_motif_within_ACR', help='Perform motif prediction using fimo tool.'
                             'Default: yes')

    ##step 02 prepare the data for the enrichment
    parser.add_argument("-s2_open_prepare_data", dest = 'open_prepare_data_enrichment', help = 'Prepare data before performing the enrichment test.'
                                                                                            'Default: yes')

    ##step 03 perform the enrichment test
    parser.add_argument("-s3_open_enrichment", dest = 'open_enrichment_test', help = 'Perform the motif enrichment test.'
                                                                                  'Default: yes')


    ##optional parameters
    ##step 01
    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument("-qval", dest='q_value_cutoff_fimo', help = 'q value cutoff for the fimo.'
                                                                    'Default: 0.00005')

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

    ##updating 040925
    ##for the chromVAR
    parser.add_argument("-acr_Tn5_fl", dest = 'acr_Tn5_file', help = 'An ACR sparse file with Tn5 integration site per peak per cell based on binary distribution.'
                                                                      'This file is generated from the step -open_motif_enrich')

    parser.add_argument("-BSgenome", dest = 'BSgenome_name', help = 'Provide a BSgenome name for the built genome of a species.')

    parser.add_argument("-jasmotifmtx", dest = 'jasmotif_matrix', help = 'Provide a jasmotif matrix. This file is obtained from the github.')

    ##the chromVAR needs the acr and meta files

    ##option 042425
    ##add an option to smooth the deviation score
    parser.add_argument("-open_smooth_dev_score", dest = 'open_smooth_motif_dev_score' ,help = 'Smooth the deviation score.'
                                                                                               'Default: yes')

    parser.add_argument("-soc_obj", dest='soc_object_fl', help='Provide an object obtained from the clustering step.')


    ##updating 050725
    ##add a function to plot the smooth motif deviation in the UMAP





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

    ##default setting
    open_predict_motif_within_ACR_final = 'no'
    open_prepare_data_enrichment_final = 'no'
    open_enrichment_test_final = 'no'
    down_sample_cell_final = 'yes'
    open_smooth_motif_dev_score_final = 'yes'


    ##motif enrichment test
    if args.open_motif_enrichment_analysis is not None:

        if args.open_motif_enrichment_analysis == 'yes':

            open_motif_enrichment_analysis_final = 'yes'

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

                #if args.genome_fai_file is None:
                #    print('Cannot find genome fai index file, please provide it')
                #    return
                #else:
                #    try:
                #        file = open(args.genome_fai_file, 'r')  ##check if the file is not the right file
                #    except IOError:
                #        print('There was an error opening the fai index file!')
                #        return

            if open_enrichment_test_final == 'yes':

                if args.celltype_cluster_col_name is None:
                    print('Please provide the name of column specifying the cell type identity in the meta file.')



        else:

            if args.open_motif_enrichment_analysis == 'no':

                open_motif_enrichment_analysis_final = 'no'

                print("Users can use \'-open_motif_enrich yes\', to open this step")
            else:
                print("Users can only use the yes|no to open or close the motif enrichment analysis.")
                return

    else:

        open_motif_enrichment_analysis_final = 'no'

        #print("Users must use \'-open_motif_enrich yes|no\', to decide if they open or close this step.")
        #return


    ##motif deviation score
    if args.open_motif_deviation_score is not None:

        if args.open_motif_deviation_score == 'yes':

            open_motif_deviation_score_final = 'yes'

            if args.acr_Tn5_file is None:
                print('Cannot find ACR Tn5 file, please provide it')
                return
            else:
                try:
                    file = open(args.acr_Tn5_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the ACR Tn5 file!')
                    return

            if args.BSgenome_name is None:
                print('please provide the BSgenome name that will be loaded from R')
                return

            if args.meta_file is None:
                print('Cannot find cell meta file, please provide it')
                return
            else:
                try:
                    file = open(args.meta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the cell meta file!')
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

            if args.jasmotif_matrix is None:
                print('Cannot find jasmotif matrix file, please provide it')
                return
            else:
                try:
                    file = open(args.jasmotif_matrix, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the jasmotif matrix file!')
                    return

            ##updating 042425
            if args.open_smooth_motif_dev_score is None:
                open_smooth_motif_dev_score_final = 'yes'
            else:
                if args.open_smooth_motif_dev_score == 'yes':
                    open_smooth_motif_dev_score_final = 'yes'
                else:
                    if args.open_smooth_motif_dev_score == 'no':
                        open_smooth_motif_dev_score_final = 'no'
                    else:
                        print("Please use \'-open_smooth_dev_score no\', to close this step")
                        return

            if open_smooth_motif_dev_score_final == 'yes':

                if args.soc_object_fl is None:
                    print('Cannot find soc object file from clustering step, please provide it')
                    return
                else:
                    try:
                        file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the soc object file from clustering step!')
                        return

                if args.celltype_cluster_col_name is None:
                    print('Please provide the name of column specifying the cell type identity in the meta file.')

                    return


        else:

            if args.open_motif_deviation_score == 'no':

                open_motif_deviation_score_final = 'no'

                print("Users can use \'-open_motif_DevScore yes\', to open this step")
            else:
                print("Users can only use the yes|no to open or close the motif deviation score generation.")
                return

    else:

        open_motif_deviation_score_final = 'no'
        #print("Users can use \'-open_motif_DevScore yes\', to open this step")


    ##upating 042425





    ##build the BS genome
    if args.open_build_BSgenome_lib is not None:

        if args.open_build_BSgenome_lib == 'yes':
            open_build_BSgenome_lib_final = 'yes'

            if args.genome_fasta_file is None:
                print('Cannot find genome fasta file, please provide it')
                return
            else:
                try:
                    file = open(args.genome_fasta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the genome fasta file!')
                    return

            if args.BSgenome_name is None:
                print('please provide the BSgenome name that will be loaded from R')
                return

        else:
            if args.open_build_BSgenome_lib == 'no':
                open_build_BSgenome_lib_final = 'no'

            else:
                print("Users can only use the yes|no to open or close the BS genome build.")
                return

    else:
        open_build_BSgenome_lib_final = 'no'


    ##other parameters
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

    ###################################################
    ##If users choose to open the motif enrichment test
    if open_motif_enrichment_analysis_final == 'yes':

        step01_motif_enrichment_analysis_dir = output_dir + '/step01_motif_enrichment_analysis_dir'
        if not os.path.exists(step01_motif_enrichment_analysis_dir):
            os.makedirs(step01_motif_enrichment_analysis_dir)

        if open_predict_motif_within_ACR_final == 'yes':

            print ('Predict motif within ACRs')

            open_predict_motif_within_ACR_final_dir = step01_motif_enrichment_analysis_dir + '/open_predict_motif_within_ACR_final_dir'
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

            open_prepare_data_enrichment_final_dir = step01_motif_enrichment_analysis_dir + '/open_prepare_data_enrichment_final_dir'
            if not os.path.exists(open_prepare_data_enrichment_final_dir):
                os.makedirs(open_prepare_data_enrichment_final_dir)

            ###############
            input_peak_fl = args.acr_file
            input_tn5_bed_fl = args.tn5_file
            #input_genome_fai_fl = args.genome_fai_file
            input_fastSparsetn5_pl = input_required_scripts_dir + '/utils_motif_analysis/fastSparse.tn5.pl'
            prefix = sample_prefix_final

            s2_prepare_data.prepare_peak_tn5_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
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
            ipt_R_script = input_required_scripts_dir + '/utils_motif_analysis/s2_prepare_regression_enrich.R'
            input_motif_acr_fl = step01_motif_enrichment_analysis_dir + '/open_predict_motif_within_ACR_final_dir' + '/opt_final_ACR_motif.bed'
            input_motif_flt_cutoff = q_value_cutoff_fimo_final

            s2_prepare_data.prepare_motif_acr_data (ipt_R_script,input_motif_acr_fl,input_motif_flt_cutoff,input_peak_cell_sparse_fl,prefix,open_prepare_data_enrichment_final_dir)

        if open_enrichment_test_final == 'yes':

            print ('Perform motif enrichment test')

            open_enrichment_test_final_dir = step01_motif_enrichment_analysis_dir + '/open_enrichment_test_final_dir'
            if not os.path.exists(open_enrichment_test_final_dir):
                os.makedirs(open_enrichment_test_final_dir)

            ipt_R_script = input_required_scripts_dir + '/utils_motif_analysis/s3_motif_enrichment_test.R'

            if down_sample_cell_final == 'yes':
                ipt_final_meta_fl = step01_motif_enrichment_analysis_dir + '/open_prepare_data_enrichment_final_dir' + '/opt_downsample_meta.txt'
            else:
                prefix_meta = 'addACRnum'
                ipt_final_meta_fl = step01_motif_enrichment_analysis_dir + '/open_prepare_data_enrichment_final_dir' + '/opt_' + prefix_meta + '.txt'

            ipt_target_cluster_nm = args.celltype_cluster_col_name

            cmd = 'Rscript ' + ipt_R_script + \
                  ' ' + step01_motif_enrichment_analysis_dir + '/open_prepare_data_enrichment_final_dir' + \
                  ' ' + ipt_final_meta_fl + \
                  ' ' + open_enrichment_test_final_dir + \
                  ' ' + ipt_target_cluster_nm + \
                  ' ' + core_number_run
            print(cmd)
            subprocess.call(cmd, shell=True)


    #######################################################################################
    ##If users choose to open the ChromVAR to perform the motif deviation score calculation
    if open_motif_deviation_score_final == 'yes':

        step02_motif_deviation_score_generation_dir = output_dir + '/step02_motif_deviation_score_generation_dir'
        if not os.path.exists(step02_motif_deviation_score_generation_dir):
            os.makedirs(step02_motif_deviation_score_generation_dir)


        ipt_R_script = input_required_scripts_dir + '/utils_motif_analysis/chromVAR_analysis.R'
        print(ipt_R_script)
        ipt_peak_tn5_fl = args.acr_Tn5_file
        ipt_acr_fl = args.acr_file
        ipt_meta_fl = args.meta_file

        BSgenome_name_str = args.BSgenome_name

        ##for the BSgenome
        store_final_parameter_line_list = []
        store_final_parameter_line_list.append('library(' + BSgenome_name_str + ')')
        ##we will firstly build up the parameter setting file
        with open(step02_motif_deviation_score_generation_dir + '/load_BSgenome.config', 'w+') as opt:
            for eachline in store_final_parameter_line_list:
                opt.write(eachline + '\n')
        ipt_load_BSgenome_config_fl = step02_motif_deviation_score_generation_dir + '/load_BSgenome.config'

        ##for the GCBias
        store_final_parameter_line_list = []
        store_final_parameter_line_list.append('fragment_counts <- addGCBias(fragment_counts, genome=' + BSgenome_name_str + ')')
        with open(step02_motif_deviation_score_generation_dir + '/load_GCBias.config', 'w+') as opt:
            for eachline in store_final_parameter_line_list:
                opt.write(eachline + '\n')
        ipt_load_GCBias_config_fl = step02_motif_deviation_score_generation_dir + '/load_GCBias.config'

        ##for the matchMotifs
        store_final_parameter_line_list = []
        store_final_parameter_line_list.append(
            'motif <- matchMotifs(jaspmotifs, filtered_counts, genome = ' + BSgenome_name_str + ')')
        with open(step02_motif_deviation_score_generation_dir + '/load_matchMotifs.config', 'w+') as opt:
            for eachline in store_final_parameter_line_list:
                opt.write(eachline + '\n')
        ipt_load_matchmotif_config_fl = step02_motif_deviation_score_generation_dir + '/load_matchMotifs.config'

        ipt_jasmotif_mtx_fl = args.jasmotif_matrix

        cmd = 'Rscript ' + ipt_R_script + \
              ' ' + core_number_run + \
              ' ' + ipt_peak_tn5_fl + \
              ' ' + ipt_acr_fl + \
              ' ' + ipt_meta_fl + \
              ' ' + ipt_jasmotif_mtx_fl + \
              ' ' + step02_motif_deviation_score_generation_dir + \
              ' ' + ipt_load_BSgenome_config_fl + \
              ' ' + ipt_load_GCBias_config_fl + \
              ' ' + ipt_load_matchmotif_config_fl
        print(cmd)
        subprocess.call(cmd,shell=True)


        ##updating 042425
        if open_smooth_motif_dev_score_final == 'yes':

            ipt_svd_obj_fl = args.soc_object_fl

            ipt_motif_deviation_score_fl = step02_motif_deviation_score_generation_dir + '/' + 'motif.deviations.txt'

            print('- open conduct smooth for the motif deviation score')

            motif_smooth_script = input_required_scripts_dir + '/utils_motif_analysis/smooth_motif_deviation.R'

            ipt_target_cluster_nm = args.celltype_cluster_col_name

            cmd = 'Rscript ' + motif_smooth_script + \
                  ' ' + ipt_meta_fl + \
                  ' ' + ipt_motif_deviation_score_fl + \
                  ' ' + ipt_svd_obj_fl + \
                  ' ' + sample_prefix_final + \
                  ' ' + ipt_target_cluster_nm + \
                  ' ' + step02_motif_deviation_score_generation_dir
            print(cmd)
            subprocess.call(cmd, shell=True)



    ######################
    ##Optional: help users to build the BSgenome
    if open_build_BSgenome_lib_final == 'yes':

        step00_build_BSgenome_dir = output_dir + '/step00_build_BSgenome_dir'
        if not os.path.exists(step00_build_BSgenome_dir):
            os.makedirs(step00_build_BSgenome_dir)

        input_genome_fasta_fl = args.genome_fasta_file

        ##build the BSgenome
        src_seqdir = step00_build_BSgenome_dir + '/src_seqdir'
        if not os.path.exists(src_seqdir):
            os.makedirs(src_seqdir)
        store_seq_dic = {}

        for seq_record in SeqIO.parse(input_genome_fasta_fl, 'fasta'):
            store_seq_dic[seq_record.id] = str(seq_record.seq)

        store_chr_list = []
        for eachid in store_seq_dic:

            ##add the fa or fasta at the same time
            opt_fl_nm = eachid + '.fa'
            with open(src_seqdir + '/' + opt_fl_nm, 'w+') as opt:
                opt.write('>' + eachid + '\n' + store_seq_dic[eachid])

            opt_fl_nm = eachid + '.fasta'
            with open(src_seqdir + '/' + opt_fl_nm, 'w+') as opt:
                opt.write('>' + eachid + '\n' + store_seq_dic[eachid])

            store_chr_list.append(eachid)



        abs_src_seqdir_path = os.path.abspath(src_seqdir)

        store_chr_list_new = []
        for eachchr in store_chr_list:
            new_chr = '\'' +  eachchr + '\''
            store_chr_list_new.append(new_chr)
        store_chr_str_new = ','.join(store_chr_list_new)

        #sep_genome(input_genome_fl, input_output_dir)
        store_final_line_list = []
        BSgenome_ID = args.BSgenome_name
        final_line = 'Package: ' + BSgenome_ID + '\n' + \
                     'Title: Full genome of ' + BSgenome_ID + '\n' + \
                     'Description: Full genome of ' + BSgenome_ID + '\n' + \
                     'Version: 1.0.0' + '\n' + \
                     'organism: ' + BSgenome_ID + '\n' + \
                     'Author: Genome' + '\n' + \
                     'common_name: Genome' + '\n' + \
                     'provider: Genome' + '\n' + \
                     'provider_version: Acyr_1.0' + '\n' + \
                     'release_date: 0.0.0' + '\n' + \
                     'release_name: ' + BSgenome_ID + '\n' + \
                     'BSgenomeObjname: Genome' + '\n' + \
                     'source_url: xxxxx' + '\n' + \
                     'organism_biocview: AnnotationData, BSgenome' + '\n' + \
                     'seqnames: c(' + store_chr_str_new + ')' + '\n' + \
                     'SrcDataFiles: Genome' + '\n' + \
                     'seqs_srcdir: ' + abs_src_seqdir_path
        store_final_line_list.append(final_line)

        with open (step00_build_BSgenome_dir + '/seed_file.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        ##it does not work for the forge script so we will not help build the forage

        ##use the R to build the seed.dcf
        #ipt_build_BSgenome_R_script = input_required_scripts_dir + '/utils_motif_analysis/build_BSgenome/forge_BSgenome.R'

        #cmd = 'Rscript ' + ipt_build_BSgenome_R_script + \
        #      ' ' + step00_build_BSgenome_dir + '/seed_file.txt' + \
        #      ' ' + step00_build_BSgenome_dir + \
        #      ' ' + BSgenome_ID
        #print(cmd)
        #subprocess.call(cmd,shell=True)

        ##use the linux to run
        #cmd = 'R CMD build ' + BSgenome_ID
        #print(cmd)
        #subprocess.call(cmd,shell=True)

        #cmd = 'R CMD check ' + BSgenome_ID + '_1.0.0.tar.gz'
        #print(cmd)
        #subprocess.call(cmd,shell=True)

        #cmd = 'R CMD INSTALL ' + BSgenome_ID + '_1.0.0.tar.gz'
        #print(cmd)
        #subprocess.call(cmd,shell=True)






if __name__ == "__main__":
    main()


















