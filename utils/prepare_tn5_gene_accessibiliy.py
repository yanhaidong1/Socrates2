#!/usr/bin/env python

import argparse
import glob
import sys
import os
import re
import subprocess


##this script is to process data bam file from cellranger
##calculate tn5 data and gene accessibility
from input_other_required_scripts_dir.utils_prepare_tn5_gene_accessibility.s1_open_process_BAM_final import pipeline_process_bam as s1_procbam
from input_other_required_scripts_dir.utils_prepare_tn5_gene_accessibility.s2_open_prepare_gene_tn5_final import pipeline_prepare_gene_tn5 as s2_geneTn5
from input_other_required_scripts_dir.utils_prepare_tn5_gene_accessibility.s3_open_gene_accessibility_final import pipeline_calculate_gene_accessibility as s3_geneAcc
from input_other_required_scripts_dir.utils_prepare_tn5_gene_accessibility.s4_open_smooth_gene_accessibility_final import pipeline_smooth_gene_accessibility as s4_geneSmooth


def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##step 01
    parser.add_argument("-BAM_fl", dest='bam_file', help="Provide bam file obtained from cell ranger")

    ##step 02
    ##updating 010425
    parser.add_argument("-soc_obj", dest = 'soc_object_fl', help ='Provide an object obtained from the clustering step.')


    parser.add_argument("-Ggff_fl", dest='gene_gff_file', help='Provide the gene gff file used to extract exon region when permutating peaks.')

    #parser.add_argument("-Gmfai_fl", dest='genome_fai_file', help='Provide a fai index file of reference genome.')

    parser.add_argument("-Genome_fl", dest = 'genome_fasta_file', help='Povide a genome fasta file.')

    ##output from the step 01
    parser.add_argument("-tn5_fl", dest = 'tn5_file', help = 'Provide a tn5 insertion file obtained from the -open_procBAM step.')


    ##step 03
    ##Users could directly use the soc object generated from the last step
    #parser.add_argument("-meta_fl", dest='meta_file', help='Provide a meta file aftering the clustering')

    ##output from the step 02
    #parser.add_argument("-gene_tn5_fl", dest='gene_tn5_file', help = 'Provide a gene tn5 file obtained from the -open_geneTn5 step.')

    #parser.add_argument("-gene_bed_fl", dest = 'gene_bed_file', help = 'Provide a gene bed file obtained from the -open_geneTn5 step.')

    ##step 04
    #parser.add_argument("-svd_fl", dest = 'svd_file', help='Provide a svd file aftering the clusteirng')

    parser.add_argument("-clustnm", dest = 'target_cluster_nm', help = 'Provide a target cluster name in the meta file to smooth the accessibility.'
                                                                       'Default: LouvainClusters')

    ##output from the step03
    #parser.add_argument("-gene_acc_fl", dest = 'gene_accessibility_file', help = 'Provide a gene accessibility file obtained from the -open_geneAcc step.')



    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')



    ##Optional parameters
    parser.add_argument("-open_procBAM", dest='s1_open_process_BAM',
                        help='Run the step 01 to process the BAM file to obtain the tn5.'
                             'Default: no')

    parser.add_argument("-open_geneTn5", dest='s2_open_prepare_gene_tn5', help='Run the step 02 to prepare the gene tn5 file.'
                                                                        'Default: no')

    parser.add_argument("-open_geneAcc", dest='s3_open_gene_accessibility',
                        help='Run the step 03 to calculate the gene accessibility.'
                             'Default: no')

    parser.add_argument("-open_smoothGeneAcc", dest= 's4_open_smooth_gene_accessibility',
                        help="Run the step 04 to smooth the gene accessibility."
                             "Default: no")


    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument("-mpq", dest='mapping_quality_score',help ='Provide a mapping quality score.'
                                                                   'Default: 10')

    parser.add_argument("-rmS1temp",dest='remove_step_s1_temp', help='Remove the temporary file in step01 to save the space.'
                                                                          'Default: yes')

    parser.add_argument("-category",dest='category_gff', help = 'Set gene or transcript to be calculated.'
                                                        'Default: gene')

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

    ##step01
    if args.s1_open_process_BAM is None:
        s1_open_process_BAM_final = 'no'
    else:
        if args.s1_open_process_BAM == 'yes':
            s1_open_process_BAM_final = 'yes'

            ##check the bam file
            if args.bam_file is None:
                print('Cannot find bam file, please provide it')
                return
            else:
                try:
                    file = open(args.bam_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the tn5 file!')
                    return

        else:
            s1_open_process_BAM_final = 'no'
            print('Users choose to close the BAM processing, please use \'-open_procBAM yes\' to open this step')




    ##step02
    if args.s2_open_prepare_gene_tn5 is None:
        s2_open_prepare_gene_tn5_final = 'no'
    else:
        if args.s2_open_prepare_gene_tn5 == 'yes':

            ###########################################################
            ##this object should be provided in all the following steps
            if args.soc_object_fl is None:
                print('Cannot find the soc object file, please provide it')
                return
            else:
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the soc object file!')
                    return

            s2_open_prepare_gene_tn5_final = 'yes'

            ##check the gff file
            if args.gene_gff_file is None:
                print('Cannot find gene gff file, please provide it')
                return
            else:
                try:
                    file = open(args.gene_gff_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the gene gff file!')
                    return

            if args.genome_fasta_file is None:
                print('Cannot find genome file, please provide it')
                return
            else:
                try:
                    file = open(args.genome_fasta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the genome file!')
                    return

            if args.tn5_file is None:
                print('Cannot find the tn5 insertion file, please provide it')
                return
            else:
                try:
                    file = open(args.tn5_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the tn5 insertion file!')
                    return




        else:
            s2_open_prepare_gene_tn5_final = 'no'
            print('Users choose to close the preparation of gene tn5 file, please use \'-open_geneTn5 yes\' to open this step')


    ##step 03
    if args.s3_open_gene_accessibility is None:
        s3_open_gene_accessibility_final = 'no'
    else:
        if args.s3_open_gene_accessibility == 'yes':
            s3_open_gene_accessibility_final = 'yes'
        else:
            s3_open_gene_accessibility_final = 'no'
            print('Users choose to close the gene accessibility calculation, please use \'-open_geneAcc yes\' to open this step')

    if args.s4_open_smooth_gene_accessibility is None:
        s4_open_smooth_gene_accessibility_final = 'no'
    else:
        if args.s4_open_smooth_gene_accessibility == 'yes':
            s4_open_smooth_gene_accessibility_final = 'yes'
        else:
            s4_open_smooth_gene_accessibility_final = 'no'
            print('Users choose to close the gene accessiblity smoothing, please use \'-smooth_geneAcc yes\' to open this step')



    if args.core_number is not None:
        core_number_run = args.core_number
    else:
        core_number_run = '1'

    if args.mapping_quality_score is not None:
        mapping_quality_val = args.mapping_quality_score
    else:
        mapping_quality_val = '10'

    if args.remove_step_s1_temp is not None:
        remove_step_s1_temp_final = args.remove_step_s1_temp
    else:
        remove_step_s1_temp_final = 'yes'

    if args.category_gff is not None:
        final_category_gff = args.category_gff
    else:
        final_category_gff = 'gene'

    if args.target_cluster_nm is not None:
        target_cluster_nm_final = args.target_cluster_nm
    else:
        target_cluster_nm_final = 'LouvainClusters'




    if s1_open_process_BAM_final == 'yes':

        print('Users choose to run the step 01 to prepare process the bam file from cell ranger')

        s1_open_process_BAM_final_dir = output_dir + '/open_process_BAM_final_dir'
        if not os.path.exists(s1_open_process_BAM_final_dir):
            os.makedirs(s1_open_process_BAM_final_dir)

        input_bam_fl = args.bam_file

        s1_procbam.process_bam (input_required_scripts_dir,input_bam_fl,s1_open_process_BAM_final_dir,
                                core_number_run,mapping_quality_val,remove_step_s1_temp_final)

    ##updating 010425 we will get the input prefix of the soc obj fl
    input_soc_obj_fl = args.soc_object_fl
    if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
        mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
        input_prefix = mt.group(1)
    else:
        if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
            mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
            input_prefix = mt.group(1)
        else:
            print('Please use *.atac.soc.rds file without changing the file name')
            return

    if s2_open_prepare_gene_tn5_final == 'yes':

        print('Users choose to run the step 02 to prepare the gene tn5 file')

        s2_open_prepare_gene_tn5_final_dir = output_dir + '/open_prepare_gene_tn5_final_dir'
        if not os.path.exists(s2_open_prepare_gene_tn5_final_dir):
            os.makedirs(s2_open_prepare_gene_tn5_final_dir)

        input_tn5_bed_fl = args.tn5_file
        #input_tn5_bed_fl = output_dir + '/open_process_BAM_final_dir/' + '/opt_tn5_mq' + mapping_quality_val + '.bed'
        input_gene_gff_fl = args.gene_gff_file
        #input_genome_fai_fl = args.genome_fai_file
        input_genome_fl = args.genome_fasta_file

        ##build the fai index
        ##updating 010425
        cmd = 'samtools faidx ' + input_genome_fl
        print(cmd)

        input_genome_fai_fl = input_genome_fl + '.fai'


        ##updating 010425
        s2_geneTn5.prepare_gene_tn5 (input_required_scripts_dir,input_tn5_bed_fl,s2_open_prepare_gene_tn5_final_dir,
                      input_gene_gff_fl,input_genome_fai_fl,
                      input_soc_obj_fl,
                      input_prefix,
                      final_category_gff)

    if s3_open_gene_accessibility_final == 'yes':

        print ('Users choose to calculate the gene accessibility')

        s3_open_gene_accessibility_final_dir = output_dir + '/open_gene_accessibility_final_dir'
        if not os.path.exists(s3_open_gene_accessibility_final_dir):
            os.makedirs(s3_open_gene_accessibility_final_dir)

        ##updating 010425

        #input_gene_sparse_fl = args.gene_tn5_file
        #input_gene_bed_fl = args.gene_bed_file

        #input_gene_sparse_fl = output_dir + '/open_prepare_gene_tn5_final_dir/opt_gene_chnm.sparse'
        #input_gene_bed_fl = output_dir + '/open_prepare_gene_tn5_final_dir/opt_genes_500bpTSS_sorted.bed'
        #input_meta_fl = args.meta_file

        s3_geneAcc.calculate_gene_body_acc(input_required_scripts_dir,
                                           input_soc_obj_fl,
                                           input_prefix,
                                           s3_open_gene_accessibility_final_dir)

    if s4_open_smooth_gene_accessibility_final == 'yes':

        print ('Users choose to smooth the gene accessibility')

        s4_open_smooth_gene_accessibility_final_dir = output_dir + '/open_smooth_gene_accessibility_final_dir'
        if not os.path.exists(s4_open_smooth_gene_accessibility_final_dir):
            os.makedirs(s4_open_smooth_gene_accessibility_final_dir)

        #input_meta_fl = args.meta_file
        #input_svd_fl = args.svd_file

        #input_gene_accessibility_mtx_rds_fl = args.gene_accessibility_file

        #input_gene_accessibility_mtx_rds_fl = output_dir + '/open_gene_accessibility_final_dir/opt_gene_body_accessibility_mtx.rd'

        s4_geneSmooth.smooth_gene_accessibility (input_required_scripts_dir,
                                                 input_soc_obj_fl,
                                                 s4_open_smooth_gene_accessibility_final_dir,core_number_run,
                                                 target_cluster_nm_final,input_prefix)

        ##remove the temp files
        cmd = 'rm ' + s4_open_smooth_gene_accessibility_final_dir + '/*sparse' + ' ' + s4_open_smooth_gene_accessibility_final_dir + '/*.txt' + ' ' + \
              s4_open_smooth_gene_accessibility_final_dir + '/all.normalizedActivity_mtx.rds'
        print(cmd)
        subprocess.call(cmd,shell=True)




if __name__ == "__main__":
    main()









