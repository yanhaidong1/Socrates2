#!/usr/bin/env python

##updating 073125 we will perform the summary of the species comparison
##updating 073025 we will perform the concise version


import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats


from s3_input_required_scripts_dir import s2_subfunctions
from s3_input_required_scripts_dir import s1_subfunctions
from s3_input_required_scripts_dir import s0_subfunctions


##########################
#input_required_scripts_dir = sys.argv[1]
##input_required_scripts_concise

########
##step00
input_spe1_syntenic_genes_all_os_ACRs_fl = sys.argv[1]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir_072925/store_results_for_next_step_dir/ipt_Pg.syntenic_genes.ACR.bed

input_all_syntenic_regions_os_acrs_dir = sys.argv[2]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir_072925/store_results_for_next_step_dir/collect_all_syntenic_regions_acrs_dir

input_all_celltype_acr_dir = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_ct_acr_Pablo_110723/final_collect_celltype_acr_dir_111623

##this is also used in the step01
input_all_spe_gene_gff_dir = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_gene_gff_111623

########
##step01
input_rice_to_allspe_blast_dir = sys.argv[5]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir_072925/store_results_for_next_step_dir/collect_output_dir_Pg_to_allspe_blast_dir

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/collect_rice_to_allspe_blast_dir_112023

input_rice_to_allspe_blast_record_riceID_dir = sys.argv[6]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir_072925/store_results_for_next_step_dir/collect_output_dir_Pg_to_allspe_blast_record_PgID_dir

input_all_syntenic_genes_all_os_ACRs_dir = sys.argv[7]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir_072925/store_results_for_next_step_dir/collect_all_other_spe_syntenic_genes_all_os_ACRs_dir

input_configure_fl = sys.argv[8]

input_output_dir = sys.argv[9]


########
##step00
##we will make a summary of how many syntenic regions

##1 How many rice syntenic regions, and what's the size of the syntenic block to the whole genome
##2 what's the distribution of distance per syntenic regions
##3 How many rice ACR in the syntenic regions (Here we would find the not shared ACRs in the rice)
##4 add these not shared ACRs to the step01
##5 step01 need to add the broad or restricted ones

def step00_basic_characters (input_output_dir,input_spe1_syntenic_genes_all_os_ACRs_fl,
                             input_all_syntenic_regions_os_acrs_dir,
                             input_all_celltype_acr_dir,
                             s0_s1_open_check_syntenic_block,s0_s1_target_spe1,s0_s1_target_spe2_str):

    step00_basic_characters_dir = input_output_dir + '/step00_basic_characters_dir'
    if not os.path.exists(step00_basic_characters_dir):
        os.makedirs(step00_basic_characters_dir)

    spe1_prefix = s0_s1_target_spe1
    input_spe1_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

    if s0_s1_open_check_syntenic_block == 'yes':

        s0_s1_open_check_syntenic_block_dir = step00_basic_characters_dir + '/s0_s1_open_check_syntenic_block_dir'
        if not os.path.exists(s0_s1_open_check_syntenic_block_dir):
            os.makedirs(s0_s1_open_check_syntenic_block_dir)

        s0_s1_target_spe2_list = s0_s1_target_spe2_str.split(',')

        ##updating 112023
        ##we will use the loop other than run several times
        for eachspe2 in s0_s1_target_spe2_list:


            ##for the rice to maize
            store_rice_to_maize_dir = s0_s1_open_check_syntenic_block_dir + '/store_' + s0_s1_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)


            spe2_prefix = eachspe2

            input_maize_syntenic_regions_os_acrs_fl = input_all_syntenic_regions_os_acrs_dir + '/' + spe2_prefix + '.syntenic_regions.os_acrs.bed'

            #input_rice_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe1_prefix + '_H3K27me3_cate.txt'

            s0_subfunctions.subfunction_check_syntenic_blocks(input_spe1_syntenic_genes_all_os_ACRs_fl, input_maize_syntenic_regions_os_acrs_fl, input_spe1_acr_fl,
                                              store_rice_to_maize_dir,
                                              spe1_prefix, spe2_prefix)

        ##for the rice to build the all the syntenic region file
        spe1_prefix = s0_s1_target_spe1
        s0_subfunctions.subfunction_build_syntenic_bed(input_spe1_syntenic_genes_all_os_ACRs_fl, input_spe1_acr_fl,s0_s1_open_check_syntenic_block_dir, spe1_prefix)





########
##step01
def step01_species_compare_add_cate (input_rice_to_allspe_blast_dir,
                                     input_all_celltype_acr_dir,
                                     input_all_spe_gene_gff_dir,
                                     input_rice_to_allspe_blast_record_riceID_dir,
                                     input_spe1_syntenic_genes_all_os_ACRs_fl,
                                     input_all_syntenic_genes_all_os_ACRs_dir,
                                     input_output_dir,s0_s1_target_spe1,s0_s1_target_spe2_str):

    step01_species_compare_add_cate_dir = input_output_dir + '/step01_species_compare_add_cate_dir'
    if not os.path.exists(step01_species_compare_add_cate_dir):
        os.makedirs(step01_species_compare_add_cate_dir)

    s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

    for eachspe2 in s0_s1_target_spe2_str_list:

        #######################
        ##for the rice to maize
        store_rice_to_maize_dir = step01_species_compare_add_cate_dir + '/store_' + s0_s1_target_spe1 + '_to_' + eachspe2 + '_dir'
        if not os.path.exists(store_rice_to_maize_dir):
            os.makedirs(store_rice_to_maize_dir)

        spe1_prefix = s0_s1_target_spe1
        spe2_prefix = eachspe2
        spe1_spe2_syntenic_acr_fl =  input_output_dir + \
                                     '/step00_basic_characters_dir/s0_s1_open_check_syntenic_block_dir' + '/store_' + spe1_prefix + '_to_' + \
                                     spe2_prefix + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_contain_ACRs.txt'

        input_rice_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

        if spe1_prefix == 'rice':
            input_rice_to_maize_blast_fl = input_rice_to_allspe_blast_dir + '/' + 'Os_ACRs.vs.' + spe2_prefix + '_blast.bed'
            input_rice_to_maize_blast_record_riceID_fl = input_rice_to_allspe_blast_record_riceID_dir + '/Os_ACRs.vs.' + spe2_prefix + '_blast_record_riceID.bed'

        else:
            ##updating 041724
            input_rice_to_maize_blast_fl = input_rice_to_allspe_blast_dir + '/' + spe1_prefix + '_ACRs.vs.' + spe2_prefix + '_blast.bed'
            ##updating 041724
            input_rice_to_maize_blast_record_riceID_fl = input_rice_to_allspe_blast_record_riceID_dir + '/' + spe1_prefix + '_ACRs.vs.' + spe2_prefix + '_blast_record_' + spe1_prefix + 'ID.bed'


        ## we do not need the H3K27me3
        ##check H3K27me3 if exist
        #input_rice_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe1_prefix + '_H3K27me3_cate.txt'
        #input_maize_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe2_prefix + '_H3K27me3_cate.txt'

        input_spe2_acr_fl = input_all_celltype_acr_dir + '/' + spe2_prefix + '_acr_celltype.txt'

        print('spe2 is ' + spe2_prefix)

        s1_subfunctions.subfunction_iden_correspond(input_rice_to_maize_blast_fl, input_rice_to_maize_blast_record_riceID_fl,
                                    input_rice_acr_fl, spe1_spe2_syntenic_acr_fl, store_rice_to_maize_dir, input_spe2_acr_fl,spe1_prefix, spe2_prefix)



        input_spe1_gff_fl = input_all_spe_gene_gff_dir + '/' + spe1_prefix + '_gene.gff'
        input_spe2_gff_fl = input_all_spe_gene_gff_dir + '/' + spe2_prefix + '_gene.gff'
        input_summary_fl = store_rice_to_maize_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_addCelltype_SynRegion.txt'

        input_maize_acr_fl = input_all_celltype_acr_dir + '/' + spe2_prefix + '_acr_celltype.txt'

        input_maize_syntenic_genes_all_os_ACRs_fl = input_all_syntenic_genes_all_os_ACRs_dir + '/' + spe2_prefix + '.syntenic_genes.all_' + spe2_prefix + '_ACRs.bed'


        s1_subfunctions.subfunction_check_syntenic_region_gene_pair_direction(input_spe1_gff_fl,input_spe2_gff_fl,
                                                                              input_spe1_syntenic_genes_all_os_ACRs_fl,input_maize_syntenic_genes_all_os_ACRs_fl,
                                                                              input_rice_to_maize_blast_fl,
                                                                              input_summary_fl, store_rice_to_maize_dir,input_maize_acr_fl,
                                                                              spe1_prefix, spe2_prefix)

        ##updating 123123
        ipt_final_summary_fl = store_rice_to_maize_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT.txt'
        #spe1_H3K27me3_fl = input_rice_H3K27me3_fl
        spe1_acr_fl = input_rice_acr_fl
        s1_subfunctions.subfunction_add_species_specific_one_under_non_syntenic_regions (ipt_final_summary_fl,spe1_acr_fl,
                                                                     spe1_prefix,spe2_prefix, store_rice_to_maize_dir)


########
##step02
def step02_summarize_overview (input_output_dir):


    ##calculate the number of syntenic regions
    step02_summarize_overview_dir = input_output_dir + '/step02_summarize_overview_dir'
    if not os.path.exists(step02_summarize_overview_dir):
        os.makedirs(step02_summarize_overview_dir)

    all_spe_pair_dir_list = glob.glob(input_output_dir + '/step01_species_compare_add_cate_dir/*')

    for eachspedir in all_spe_pair_dir_list:

        mt = re.match('.+/(.+)',eachspedir)
        dirnm = mt.group(1)
        mt = re.match('store_(.+)_to_(.+)_dir',dirnm)
        spe1nm = mt.group(1)
        spe2nm = mt.group(2)

        ipt_fl = eachspedir + '/opt_' + spe1nm + '_' + spe2nm + '_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addSpeSpec.txt'

        s2_subfunctions.subfunction_summarize_overview(ipt_fl, step02_summarize_overview_dir, spe1nm, spe2nm)


















store_target_parameter_dic = {}
with open (input_configure_fl,'r') as ipt:
    for eachline in ipt:
        eachline = eachline.strip('\n')
        if not eachline.startswith('#'):
            col = eachline.strip().split('=')
            store_target_parameter_dic[col[0]] = col[1]

step00 = store_target_parameter_dic['step00'] ##yes or no
step01 = store_target_parameter_dic['step01'] ##yes or no
step02 = store_target_parameter_dic['step02'] ##yes or no

input_core_num = store_target_parameter_dic['input_core_num']

s0_s1_open_check_syntenic_block = store_target_parameter_dic['s0_s1_open_check_syntenic_block']
s0_s1_target_spe1 = store_target_parameter_dic['s0_s1_target_spe1']
s0_s1_target_spe2_str = store_target_parameter_dic['s0_s1_target_spe2_str']
s0_s1_width_bin = store_target_parameter_dic['s0_s1_width_bin']



if step00 == 'yes':
    step00_basic_characters(input_output_dir, input_spe1_syntenic_genes_all_os_ACRs_fl,
                            input_all_syntenic_regions_os_acrs_dir,
                            input_all_celltype_acr_dir,
                            s0_s1_open_check_syntenic_block, s0_s1_target_spe1, s0_s1_target_spe2_str)

if step01 == 'yes':
    step01_species_compare_add_cate(input_rice_to_allspe_blast_dir,
                                    input_all_celltype_acr_dir,
                                    input_all_spe_gene_gff_dir,
                                    input_rice_to_allspe_blast_record_riceID_dir,
                                    input_spe1_syntenic_genes_all_os_ACRs_fl,
                                    input_all_syntenic_genes_all_os_ACRs_dir,
                                    input_output_dir, s0_s1_target_spe1, s0_s1_target_spe2_str)

if step02 == 'yes':
    step02_summarize_overview(input_output_dir)