#!/usr/bin/env python

import argparse
import glob
import sys
import os
import re
import subprocess

##updating 101825 we will add the defined an option for users to select the
##updating 031625 add the black list building
##updating 020325 add chr that not consider in the gene features
##updating 012325 add removing the black list step

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
    parser.add_argument("-BAM_fl", dest='bam_file', help="Provide bam file obtained from cell ranger.")

    parser.add_argument("-black_fl", dest='black_list_file', help = 'Provide a black list file.')


    ##################
    ##build black list pipeline update 031625
    parser.add_argument("-open_build_black", dest = 's1_open_build_black_file', help = 'Users choose to build a black list file.')

    ##option genomic black
    parser.add_argument("-genomic_black",dest = 'open_genomic_black',help = 'Users choose to use the genomic control file to build the black file.'
                                                                                 'Default:yes')

    parser.add_argument("-genomic_bed_fl", dest = 'genomic_bed_file', help = 'Users provide a genomic bed file if open the geonmic black option.')

    parser.add_argument("-Gmfai_fl", dest='genome_fai_file', help='Provide a fai index file of reference genome.')

    ##option mtpt black
    parser.add_argument("-organelle_black", dest = 'open_organelle_black',help = 'Users choose to build the organelle control file.'
                                                                                'Default: yes')

    #parser.add_argument("-Genome_fl", dest='genome_fasta_file', help='Povide a genome fasta file.')

    parser.add_argument("-organelle_chr_name", dest = 'organelle_chromosome_name', help = 'Provide a organelle chromosome name seperated by comma.'
                                                                        )

    ##option repeat black
    parser.add_argument("-repeat_black",dest = 'open_repeat_black', help = 'Users choose to build the repeat black control file.'
                                                                           'Default: yes')

    parser.add_argument("-repeatmasker_fl", dest = 'repeat_masker_file', help = 'Provide a repeat masker file.')






    ##step 02
    ##updating 010425
    parser.add_argument("-soc_obj", dest = 'soc_object_fl', help ='Provide an object obtained from the clustering step.')


    parser.add_argument("-Ggff_fl", dest='gene_gff_file', help='Provide the gene gff file used to extract exon region when permutating peaks.')

    parser.add_argument("-black_chr",  dest = 'black_chromosome_string', help = 'Provide a list of chromosomes not considering in the features.'
                                                                         'Default: ChrPt,ChrMt')

    #parser.add_argument("-Gmfai_fl", dest='genome_fai_file', help='Provide a fai index file of reference genome.')

    parser.add_argument("-Genome_fl", dest = 'genome_fasta_file', help='Povide a genome fasta file.')

    ##output from the step 01
    parser.add_argument("-tn5_fl", dest = 'tn5_file', help = 'Provide a tn5 insertion file obtained from the -open_procBAM step.')


    ##updating 101825
    parser.add_argument("-promoter_ext_size", dest = 'promoter_extended_size', help ='Users will define a range of using exact promoter extended size for calculating the gene chromatin accessibility. Default is 500 bp'
                                                                                     'Default: 500')



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

    ##updating 010526
    parser.add_argument("-use_picard", dest = 'use_picard_tool', help ='If the users will use the picrad to remove the PCR duplicate.'
                                                                       'Default: yes')




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
                    print('There was an error opening the BAM file!')
                    return

        else:
            s1_open_process_BAM_final = 'no'
            print('Users choose to close the BAM processing, please use \'-open_procBAM yes\' to open this step')

    ##updaitng 012325
    #if args.black_list_file is None:
    #    print('Cannot find black list file, if users want to remove the reads under the black region, please provide it')
    #else:
    #    try:
    #        file = open(args.black_list_file, 'r')  ##check if the file is not the right file
    #    except IOError:
    #        print('There was an error opening the black list file!')
    #        return

    ##updating 031625
    ##make the black list file
    open_genomic_black_final = 'yes'
    open_organelle_black_final = 'yes'
    open_repeat_black_final = 'yes'
    if args.s1_open_build_black_file is None:
        s1_open_build_black_file_final = 'no'
    else:
        if args.s1_open_build_black_file == 'yes':

            s1_open_build_black_file_final = 'yes'

            ##check the genomic black
            if args.open_genomic_black is None:
                open_genomic_black_final = 'yes'
            else:
                if args.open_genomic_black == 'yes':
                    open_genomic_black_final = 'yes'
                else:
                    open_genomic_black_final = 'no'

            if open_genomic_black_final == 'yes':

                ##check the genomic bed file
                if args.genomic_bed_file is None:
                    print('Cannot find genomic bed file, please provide it')
                    return
                else:
                    try:
                        file = open(args.genomic_bed_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the genomic bed file!')
                        return

                ##check the genome_fai_file
                if args.genome_fai_file is None:
                    print('Cannot find genome fai file, please provide it')
                    return
                else:
                    try:
                        file = open(args.genome_fai_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the genome fai file!')
                        return

            ##check the organelle black
            if args.open_organelle_black is None:
                open_organelle_black_final = 'yes'
            else:
                if args.open_organelle_black == 'yes':
                    open_organelle_black_final = 'yes'
                else:
                    open_organelle_black_final = 'no'

            if open_organelle_black_final == 'yes':

                if args.genome_fasta_file is None:
                    print('Cannot find genome fasta file, please provide it')
                    return
                else:
                    try:
                        file = open(args.genome_fasta_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the genome fasta file!')
                        return

                if args.organelle_chromosome_name is None:
                    print('Users must provide the chromosome name of organelle seperated by comma')
                    return


                    ##check the repeat black
            if args.open_repeat_black is None:
                open_repeat_black_final = 'yes'
            else:
                if args.open_repeat_black == 'yes':
                    open_repeat_black_final = 'yes'
                else:
                    open_repeat_black_final = 'no'

            if open_repeat_black_final == 'yes':

                if args.repeat_masker_file is None:
                    print('Cannot find repeat masker file, please provide it')
                    return
                else:
                    try:
                        file = open(args.repeat_masker_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the repeat masker file!')
                        return

        else:
            s1_open_build_black_file_final = 'no'

            print('Users do not build the black list file, please use \'-open_build_black yes\' to open this step')




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


    ##updating 101825
    if args.promoter_extended_size is not None:
        promoter_extended_size_final = args.promoter_extended_size
    else:
        promoter_extended_size_final = '500'


    ##updating 010526
    if args.use_picard_tool is None:
        use_picard_tool_final = 'yes'
    else:
        if args.use_picard_tool == 'yes':
            use_picard_tool_final = 'yes'
        else:
            if args.use_picard_tool == 'no':
                use_picard_tool_final = 'no'
            else:
                print('Please set -use_picard no to close the picard usage')
                return



    ##updating 031625
    if s1_open_build_black_file_final == 'yes':

        open_build_black_file_final_dir = output_dir + '/open_build_black_file_final_dir'
        if not os.path.exists(open_build_black_file_final_dir):
            os.makedirs(open_build_black_file_final_dir)

        store_final_black_file_dir = open_build_black_file_final_dir + '/store_final_black_file_dir'
        if not os.path.exists(store_final_black_file_dir):
            os.makedirs(store_final_black_file_dir)

        print('Users choose to build the black file which records regions that will be filtered out in the following analysis')

        if open_genomic_black_final == 'yes':

            print ('Users will build the black file based on the genomic bed file')

            genomic_black_dir = open_build_black_file_final_dir + '/genomic_black_dir'
            if not os.path.exists(genomic_black_dir):
                os.makedirs(genomic_black_dir)

            ##step01
            ipt_intersect_python_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/genomic_black/01_intersect_with_1kwin.py'

            ipt_genomic_fl = args.genomic_bed_file
            ipt_genome_fai_fl = args.genome_fai_file
            #ipt_nobinary_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/genomic_black/fastSparse.nonbinary.peak.pl'
            ##updating 052025
            ipt_nobinary_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/genomic_black/fastSparse.nonbinary.peak.py'

            cmd = 'python ' + ipt_intersect_python_script + \
                  ' ' + ipt_genomic_fl + \
                  ' ' + ipt_genome_fai_fl + \
                  ' ' + ipt_nobinary_script + \
                  ' ' + genomic_black_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

            ##step02
            ipt_identify_black_region_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/genomic_black/02_identify_black_region.py'
            ipt_genomic_bed_fl = genomic_black_dir + '/opt_nb_win_read.sparse'

            cmd = 'python ' + ipt_identify_black_region_script + \
                  ' ' + ipt_genomic_bed_fl + \
                  ' ' + genomic_black_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

            cmd = 'cp ' + genomic_black_dir + '/opt_genomic_black.bed ' + store_final_black_file_dir
            print(cmd)
            subprocess.call(cmd,shell=True)


        if open_organelle_black_final == 'yes':

            print('Users will build the organelle black file')

            organelle_black_dir = open_build_black_file_final_dir + '/organelle_black_dir'
            if not os.path.exists(organelle_black_dir):
                os.makedirs(organelle_black_dir)

            ##Step01
            ipt_blast_python_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/mtpt_black/01_blast_mtpt_genome.py'
            ipt_genome_fl = args.genome_fasta_file
            organelle_chromosome_name_final_str = args.organelle_chromosome_name

            cmd = 'python ' + ipt_blast_python_script + \
                  ' ' + ipt_genome_fl + \
                  ' ' + organelle_black_dir + \
                  ' ' + organelle_chromosome_name_final_str
            print(cmd)
            subprocess.call(cmd,shell=True)

            ##Step02
            ipt_generate_black_script =  input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/mtpt_black/02_generate_black.py'

            cmd = 'python ' + ipt_generate_black_script + \
                  ' ' + organelle_black_dir + '/opt_blast.txt' + \
                  ' ' + organelle_black_dir + \
                  ' ' + organelle_chromosome_name_final_str
            print(cmd)
            subprocess.call(cmd,shell=True)

            cmd = 'cp ' + organelle_black_dir + '/opt_mtpt_black.bed ' + store_final_black_file_dir
            print(cmd)
            subprocess.call(cmd,shell=True)


        if open_repeat_black_final == 'yes':

            print ('Users will build the repeat black file')

            repeat_black_dir = open_build_black_file_final_dir + '/repeat_black_dir'
            if not os.path.exists(repeat_black_dir):
                os.makedirs(repeat_black_dir)

            ##Step01
            ipt_make_black_repeatmasker_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_build_black_final/repeat_black/make_black_repeatmask.py'

            ipt_repeatmasker_fl = args.repeat_masker_file

            cmd = 'python ' + ipt_make_black_repeatmasker_script + \
                  ' ' + ipt_repeatmasker_fl + \
                  ' ' + repeat_black_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

            cmd = 'cp ' + repeat_black_dir + '/opt_repeat_black.bed ' + store_final_black_file_dir
            print(cmd)
            subprocess.call(cmd, shell=True)


        cmd = 'cat ' + store_final_black_file_dir + '/* > ' + store_final_black_file_dir + '/opt_final_black_region.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)



    if s1_open_process_BAM_final == 'yes':

        print('Users choose to run the step 01 to prepare process the bam file from cell ranger')

        s1_open_process_BAM_final_dir = output_dir + '/open_process_BAM_final_dir'
        if not os.path.exists(s1_open_process_BAM_final_dir):
            os.makedirs(s1_open_process_BAM_final_dir)

        input_bam_fl = args.bam_file

        s1_procbam.process_bam (input_required_scripts_dir,input_bam_fl,s1_open_process_BAM_final_dir,
                                core_number_run,mapping_quality_val,remove_step_s1_temp_final,use_picard_tool_final)

        ##updating 012325 update the black list file
        if args.black_list_file is not None:

            print ('Users choose to use black list file to filter out Tn5 reads within black region in genome')

            ipt_black_list_file = args.black_list_file
            input_qual_val = mapping_quality_val

            ipt_tn5_fl = s1_open_process_BAM_final_dir + '/opt_tn5_mq' + input_qual_val + '.bed'

            s1_procbam.remove_black (ipt_tn5_fl,ipt_black_list_file,s1_open_process_BAM_final_dir,remove_step_s1_temp_final)



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
        subprocess.call(cmd,shell=True)
        print(cmd)

        input_genome_fai_fl = input_genome_fl + '.fai'

        if args.black_chromosome_string is None:
            final_black_chr_str = 'ChrPt,ChrMt'
        else:
            final_black_chr_str = args.black_chromosome_string

        ##updating 010425
        s2_geneTn5.prepare_gene_tn5 (input_required_scripts_dir,input_tn5_bed_fl,s2_open_prepare_gene_tn5_final_dir,
                      input_gene_gff_fl,final_black_chr_str,input_genome_fai_fl,
                      input_soc_obj_fl,
                      input_prefix,promoter_extended_size_final,
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









