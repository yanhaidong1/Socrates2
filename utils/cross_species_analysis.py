#!/usr/bin/env python

import argparse
import glob
import sys
import subprocess
import os
import re


##need to install GENESPACE in R and orthofinder in the environment

##updating 070525 we will use the object


def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    #########
    ##step 01
    parser.add_argument("-s1_open_prepare", dest='s1_open_prepare_synhit_file',
                        help='Run the step 01 to prepare syntenic region file.'
                             'Default: yes')

    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    parser.add_argument("-bed_dir", dest = 'bed_directory',help = 'Users must provide the bed dir that includes the coordinates of genes.')

    parser.add_argument("-peptide_dir", dest = 'peptide_directory',help = 'Users must provide the peptide dir that includes the protein sequences.')

    #parser.add_argument("-path2mcscanx", dest = 'path2mcscanx_path',help = 'Users must provide a ptah to path2mcscanx.')

    parser.add_argument("-ref_spe", dest = 'reference_species', help = 'Users must specify which species need to be considered as the reference species.')


    #########
    ##step 02
    parser.add_argument('-s2_open_blast_ACR', dest = 's2_open_blast_ACR_file', help = 'Run the step 02 to blast ACR to the syntenic regions.'
                                                                                      'Default: yes.')

    parser.add_argument('-cell_type_dir', dest = 'cell_type_directory', help = 'This directory contains the cell type file for all species.')

    parser.add_argument('-genome_fa_dir', dest = 'genome_fasta_dir', help = 'This directory contains the genome fasta file for all species.')


    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    #########
    ##step 03
    parser.add_argument('-gene_gff_dir', dest='gene_gff_directory',
                        help='This directory contains the gene gff files for all species.')

    parser.add_argument('-s3_open_summary', dest='s3_open_summary_ACR_syn',
                        help='Run the step 03 to summarize ACR to the syntenic regions.'
                             'Default: yes.')





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

    #if args.celltype_cluster_col_name is None:
    #    print('Please provide the name of column specifying the cell type identity in the meta file.')
    #    return

    input_required_scripts_dir = args.required_script_dir

    if args.s1_open_prepare_synhit_file is None:
        s1_open_prepare_syn_fl_final = 'yes'
    else:
        if args.s1_open_prepare_synhit_file == 'yes':
            s1_open_prepare_syn_fl_final = 'yes'
        else:
            s1_open_prepare_syn_fl_final = 'no'
            print(
                'Users choose to close the preparation of syntenic file, please use \'-s1_open_prepare yes\' to open this step')

    if s1_open_prepare_syn_fl_final == 'yes':

        if args.bed_directory is None:
            print('Cannot find required bed dir, please provide the dir in \'-bed_dir\' !')
            return
        else:
            if not os.path.exists(args.bed_directory):
                print("Bed directory does not exist, please provide it")
                return

            elif not os.listdir(args.bed_directory):
                print("Directory is empty, please make sure this directory contains the gene coordinate files")


        if args.peptide_directory is None:
            print('Cannot find required peptide dir, please provide the dir in \'-peptide_dir\' !')
            return
        else:
            if not os.path.exists(args.peptide_directory):
                print("peptide directory does not exist, please provide it")
                return

            elif not os.listdir(args.peptide_directory):
                print("Directory is empty, please make sure this directory contains the protein fasta files")


        #if args.path2mcscanx_path is None:
        #    print ('Cannot find the path of mcscanx, please provide it')
        #    return

        if args.reference_species is None:
            print('Cannot find the reference species symbol, please provide it')
            return

    if args.s2_open_blast_ACR_file is None:
        s2_open_blast_ACR_file_final = 'yes'
    else:
        if args.s2_open_blast_ACR_file == 'yes':
            s2_open_blast_ACR_file_final = 'yes'
        else:
            s2_open_blast_ACR_file_final = 'no'
            print(
                'Users choose to close the blasting ACR, please use \'-s2_open_blast_ACR yes\' to open this step')


    if s2_open_blast_ACR_file_final == 'yes':

        if args.cell_type_directory is None:
            print('Cannot find required cell type directory, please provide the dir in \'-cell_type_dir\' !')
            return
        else:
            if not os.path.exists(args.cell_type_directory):
                print("peptide directory does not exist, please provide it")
                return

            elif not os.listdir(args.cell_type_directory):
                print("Directory is empty, please make sure this dir")

        if args.genome_fasta_dir is None:
            print('Cannot find required genome fasta directory, please provide the dir in \'-genome_fa_dir\' !')
            return
        else:
            if not os.path.exists(args.genome_fasta_dir):
                print("peptide directory does not exist, please provide it")
                return

            elif not os.listdir(args.genome_fasta_dir):
                print("Directory is empty, please make sure this dir")

    if args.s3_open_summary_ACR_syn is None:
        s3_open_summary_ACR_syn_final = 'yes'
    else:
        if args.s3_open_summary_ACR_syn == 'yes':
            s3_open_summary_ACR_syn_final = 'yes'
        else:
            s3_open_summary_ACR_syn_final = 'no'
            print(
                'Users choose to close the summarizing results, please use \'-s3_open_summary yes\' to open this step')

    if s3_open_summary_ACR_syn_final == 'yes':

        if args.gene_gff_directory is None:
            print('Cannot find required genome fasta directory, please provide the dir in \'-genome_fa_dir\' !')
            return
        else:
            if not os.path.exists(args.gene_gff_directory):
                print("peptide directory does not exist, please provide it")
                return

            elif not os.listdir(args.gene_gff_directory):
                print("Directory is empty, please make sure this dir")

        if args.cell_type_directory is None:
            print('Cannot find required cell type directory, please provide the dir in \'-cell_type_dir\' !')
            return
        else:
            if not os.path.exists(args.cell_type_directory):
                print("peptide directory does not exist, please provide it")
                return

            elif not os.listdir(args.cell_type_directory):
                print("Directory is empty, please make sure this dir")



    if s1_open_prepare_syn_fl_final == 'yes':

        ########################
        ##step 01 run gene space
        step01_run_genespace_dir = output_dir + '/step01_run_genespace_dir'
        if not os.path.exists(step01_run_genespace_dir):
            os.makedirs(step01_run_genespace_dir)

        working_dir = step01_run_genespace_dir + '/working_dir'
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        ipt_bed_dir = args.bed_directory
        ipt_peptide_dir = args.peptide_directory
        #ipt_path2mcscanx = args.path2mcscanx_path


        cmd = 'cp -r ' + ipt_bed_dir + ' ' + working_dir + '/' + 'bed'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'cp -r ' + ipt_peptide_dir + ' ' + working_dir + '/' + 'peptide'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ipt_genespace_R_script = input_required_scripts_dir + '/utils_cross_species/run_GENESPACE.R'

        cmd = 'Rscript ' + ipt_genespace_R_script + \
              ' ' + working_dir
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##run the downstream
        ipt_modi_fl_python_script = input_required_scripts_dir + '/utils_cross_species/syntenic_region_pipeline.py'

        ipt_ref_spe_symbol = args.reference_species

        ipt_spe2spe_fl_list = glob.glob(working_dir + '/results/' + ipt_ref_spe_symbol + '__v__*')

        for eachspe2spe_fl in ipt_spe2spe_fl_list:

            ipt_target_spe2spe_fl = eachspe2spe_fl
            mt = re.match('.+/(.+)',ipt_target_spe2spe_fl)
            flnm = mt.group(1)

            mt = re.match('(.+)__v__(.+)\.tsv',flnm)
            sp1 = mt.group(1)
            sp2 = mt.group(2)

            if os.path.exists(working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt.gz'):
                cmd = 'gunzip ' + working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt.gz'
                print(cmd)
                subprocess.call(cmd,shell=True)
                ipt_target_syn_fl = working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt'
            else:
                if os.path.exists(working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt'):
                    ipt_target_syn_fl = working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt'
                else:

                    if os.path.exists(working_dir + '/syntenicHits/' + sp2 + '_vs_' + sp1 + '.synHits.txt.gz'):
                        cmd = 'gunzip ' + working_dir + '/syntenicHits/' + sp2 + '_vs_' + sp1 + '.synHits.txt.gz'
                        print(cmd)
                        subprocess.call(cmd, shell=True)

                        ipt_modi_syn_fl_script = input_required_scripts_dir + '/utils_cross_species/shift_species_synhits.py'

                        cmd = 'python ' + ipt_modi_syn_fl_script + \
                              ' ' + working_dir + '/syntenicHits/' + sp2 + '_vs_' + sp1 + '.synHits.txt' + \
                              ' ' + working_dir + '/syntenicHits/' + \
                              ' ' + sp1 + '_vs_' + sp2 + '.synHits'
                        print(cmd)
                        subprocess.call(cmd,shell=True)

                        ipt_target_syn_fl = working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt'

                    else:
                        if os.path.exists(working_dir + '/syntenicHits/' + sp2 + '_vs_' + sp1 + '.synHits.txt'):

                            ipt_modi_syn_fl_script = input_required_scripts_dir + '/utils_cross_species/shift_species_synhits.py'

                            cmd = 'python ' + ipt_modi_syn_fl_script + \
                                  ' ' + working_dir + '/syntenicHits/' + sp2 + '_vs_' + sp1 + '.synHits.txt' + \
                                  ' ' + working_dir + '/syntenicHits/' + \
                                  ' ' + sp1 + '_vs_' + sp2 + '.synHits'
                            print(cmd)
                            subprocess.call(cmd, shell=True)

                            ipt_target_syn_fl = working_dir + '/syntenicHits/' + sp1 + '_vs_' + sp2 + '.synHits.txt'

                        else:
                            print('Cannot find the correct synHits file')
                            return


            cmd = 'python ' + ipt_modi_fl_python_script + \
                  ' --file1 ' + ipt_target_spe2spe_fl + \
                  ' --file2 ' + ipt_target_syn_fl
            print(cmd)
            subprocess.call(cmd,shell=True)


    if s2_open_blast_ACR_file_final == 'yes':

        #######################
        ##step 02 run blast ACR
        step02_run_ACR_blast_dir = output_dir + '/step02_run_ACR_blast_dir'
        if not os.path.exists(step02_run_ACR_blast_dir):
            os.makedirs(step02_run_ACR_blast_dir)

        ipt_pipeline_blast_ACR_script = input_required_scripts_dir + '/utils_cross_species/s2_pipeline_blasted_ACRs.py'

        ipt_syntenic_region_all_fl_list = glob.glob(output_dir + '/step01_run_genespace_dir/working_dir/syntenicHits/*final_converted.txt')

        for eachsynfl in ipt_syntenic_region_all_fl_list:

            mt = re.match('.+/(.+)',eachsynfl)
            flnm = mt.group(1)

            mt = re.match('(.+)_vs_(.+)\.synHits_isAnchor_TRUE_filtered_final_converted\.txt',flnm)
            spe1 = mt.group(1)
            spe2 = mt.group(2)

            opt_dir = step02_run_ACR_blast_dir + '/' + spe1 + '_' + spe2
            if not os.path.exists(opt_dir):
                os.makedirs(opt_dir)

            ipt_cell_type_dir = args.cell_type_directory
            ipt_genome_fasta_dir = args.genome_fasta_dir
            ipt_core_num = args.core_number

            cmd = 'python ' + ipt_pipeline_blast_ACR_script + \
                  ' ' + ipt_cell_type_dir + \
                  ' ' + eachsynfl + \
                  ' ' + ipt_genome_fasta_dir + \
                  ' ' + ipt_core_num + \
                  ' ' + opt_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

    if s3_open_summary_ACR_syn_final == 'yes':

        #####################
        ##step 03 summary ACR
        step03_run_summary_ACR_syn_dir = output_dir + '/step03_run_summary_ACR_syn_dir'
        if not os.path.exists(step03_run_summary_ACR_syn_dir):
            os.makedirs(step03_run_summary_ACR_syn_dir)

        ipt_script_fl = input_required_scripts_dir + '/utils_cross_species/s3_pipeline_summary_syntenic_ACRs.py'

        all_spe_dir_list = glob.glob(output_dir + '/step02_run_ACR_blast_dir/*')

        for eachspedir in all_spe_dir_list:

            mt = re.match('.+/(.+)',eachspedir)
            dirnm = mt.group(1)

            mt = re.match('(.+)_(.+)',dirnm)
            spe1 = mt.group(1)
            spe2 = mt.group(2)

            opt_dir = step03_run_summary_ACR_syn_dir + '/' + spe1 + '_' + spe2
            if not os.path.exists(opt_dir):
                os.makedirs(opt_dir)

            ipt_core_num = args.core_number

            ipt_spe1_syn_gene_acr_fl = eachspedir + '/store_results_for_next_step_dir/ipt_' + spe1 + '.syntenic_genes.ACR.bed'
            ipt_collect_all_syn_acrs_dir = eachspedir + '/store_results_for_next_step_dir/collect_all_syntenic_regions_acrs_dir'
            ipt_cell_type_dir = args.cell_type_directory
            ipt_gene_gff_dir = args.gene_gff_directory
            ipt_collect_spe1toallspe_blast_dir = eachspedir + '/store_results_for_next_step_dir/collect_output_dir_Pg_to_allspe_blast_dir'
            ipt_collect_spe1toallspe_blast_record_spe1ID_dir = eachspedir + '/store_results_for_next_step_dir/collect_output_dir_Pg_to_allspe_blast_record_PgID_dir'
            ipt_collect_all_other_spe_syn_gene_all_ACR_dir = eachspedir + '/store_results_for_next_step_dir/collect_all_other_spe_syntenic_genes_all_os_ACRs_dir'

            ##build a temp configure file
            store_final_line_list = []
            final_line = 'step00=yes' + '\n' + 'step01=yes' + '\n' + 'step02=yes' + '\n' + \
                         'input_core_num=' + ipt_core_num + '\n' + \
                         's0_s1_open_check_syntenic_block=yes' + '\n' + \
                         's0_s1_target_spe1=' + spe1 + '\n' + \
                         's0_s1_target_spe2_str=' + spe2 + '\n' + \
                         's0_s1_width_bin=5'
            store_final_line_list.append(final_line)

            with open (opt_dir + '/temp_config.txt','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            ipt_config_fl = opt_dir + '/temp_config.txt'

            cmd = 'python ' + ipt_script_fl + \
                  ' ' + ipt_spe1_syn_gene_acr_fl + \
                  ' ' + ipt_collect_all_syn_acrs_dir + \
                  ' ' + ipt_cell_type_dir + \
                  ' ' + ipt_gene_gff_dir + \
                  ' ' + ipt_collect_spe1toallspe_blast_dir + \
                  ' ' + ipt_collect_spe1toallspe_blast_record_spe1ID_dir + \
                  ' ' + ipt_collect_all_other_spe_syn_gene_all_ACR_dir + \
                  ' ' + ipt_config_fl + \
                  ' ' + opt_dir
            print(cmd)
            subprocess.call(cmd,shell=True)






if __name__ == "__main__":
    main()








































































