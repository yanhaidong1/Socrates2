#!/usr/bin/env python

##updating 053124 we will add the os syntenic gene file
##updating 050524 we will merge and filter the blasted ACRs

import re
import glob
import sys
import subprocess
import os

from s2_input_required_scripts_dir import s1_functions
from s2_input_required_scripts_dir import s2_functions
from s2_input_required_scripts_dir import s3_functions

input_ct_acr_dir = sys.argv[1]
##/scratch/hy17471/database_031724/section3_add_pearlmillet_101424/10_downstream_analysis_071725/05_cross_species_072125/input_all_spe_celltype_acr_dir


input_syntenic_region_fl = sys.argv[2]
##/scratch/hy17471/database_031724/section3_add_pearlmillet_101424/10_downstream_analysis_071725/05_cross_species_072125/Pg_vs_Zea_Sbi_Urf_Pmi_syntenic_regions.txt

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_revising_ms_041624/08_compare_ACR_root_maize_rice_042124/15_syntenic_corresponding_050324/ipt_reference_all_regions_dir/os_reference.all_species_comparisons.filtered.tsv
##this is all the spe record
##Here rice must be the spe 1

##step02
input_all_spe_genome_dir = sys.argv[3]
##/scratch/hy17471/database_031724/section3_add_pearlmillet_101424/10_downstream_analysis_071725/05_cross_species_072125/input_all_spe_fasta_dir

##step03
#input_ct_acr_dir = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_ct_acr_root_rice_maize_Pablo_053024/final_collect_celltype_acr_dir_053124

#input_configure_fl = sys.argv[5]

input_core_num = sys.argv[4]

input_output_dir = sys.argv[5]

##updating 082325
input_spe1 = sys.argv[6]

input_spe2 = sys.argv[7]

open_step01 = 'yes'

open_step02 = 'yes'


def step01_build_ACR_in_syntenic_region (input_ct_acr_dir,input_syntenic_region_fl,
                                         s1_target_spe1_str,
                                         s1_target_spe2_str,input_output_dir):

    step01_build_ACR_in_syntenic_region_dir = input_output_dir + '/step01_build_ACR_in_syntenic_region_dir'
    if not os.path.exists(step01_build_ACR_in_syntenic_region_dir):
        os.makedirs(step01_build_ACR_in_syntenic_region_dir)

    s1_target_spe2_str_list = s1_target_spe2_str.split(',')

    for eachspe2 in s1_target_spe2_str_list:

        print('##############')
        print('target spe1 is ' + s1_target_spe1_str)
        print('target spe2 is ' + eachspe2)
        print('##############')

        #######################
        ##for the rice to maize
        store_spe1_to_spe2_dir = step01_build_ACR_in_syntenic_region_dir + '/store_' + s1_target_spe1_str + '_to_' + eachspe2 + '_dir'
        if not os.path.exists(store_spe1_to_spe2_dir):
            os.makedirs(store_spe1_to_spe2_dir)

        temp_spe1_ACR_fl = input_ct_acr_dir + '/' + s1_target_spe1_str + '_acr_celltype.txt'
        temp_spe2_ACR_fl = input_ct_acr_dir + '/' + eachspe2 + '_acr_celltype.txt'

        store_final_line_list = []
        with open (temp_spe1_ACR_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

        with open (input_output_dir + '/temp_spe1_ACR_fl.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')
        cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/temp_spe1_ACR_fl.txt > ' + input_output_dir + '/temp_spe1_ACR_fl_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        store_final_line_list = []
        with open(temp_spe2_ACR_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

        with open(input_output_dir + '/temp_spe2_ACR_fl.txt', 'w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')
        cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/temp_spe2_ACR_fl.txt > ' + input_output_dir + '/temp_spe2_ACR_fl_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

        spe1_ACR_fl = input_output_dir + '/temp_spe1_ACR_fl_sorted.txt'
        spe2_ACR_fl = input_output_dir + '/temp_spe2_ACR_fl_sorted.txt'

        ##three col

        s1_functions.s1_subfunction_build_ACR_in_syn_region (spe1_ACR_fl,spe2_ACR_fl,input_syntenic_region_fl,
                                                             store_spe1_to_spe2_dir,
                                                             s1_target_spe1_str,eachspe2)

        s1_functions.s1_subfunction_build_Os_syntenic_genes_all_os_ACR_fl(input_syntenic_region_fl, store_spe1_to_spe2_dir,
                                                                          s1_target_spe1_str)


        ##the output will be chr + '\t' + str(min_loc) + '\t' + str(max_loc) + '\t' + \
        # eachsyntenicID + '\t' + 'ACRchr + '\t' + 'ACRst' + '\t' + 'ACRed'


def step02_conduct_blast_between_two_spe (input_output_dir,input_all_spe_genome_dir,s2_target_spe1_str,s2_target_spe2_str,input_core_num):

    ##we will set different cores to do the blast
    ##each blast we will divide it into different bins

    step02_conduct_blast_between_two_spe_dir = input_output_dir + '/step02_conduct_blast_between_two_spe_dir'
    if not os.path.exists(step02_conduct_blast_between_two_spe_dir):
        os.makedirs(step02_conduct_blast_between_two_spe_dir)

    s2_target_spe2_str_list = s2_target_spe2_str.split(',')

    for eachspe2 in s2_target_spe2_str_list:

        store_spe1_to_spe2_dir = step02_conduct_blast_between_two_spe_dir + '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir'
        if not os.path.exists(store_spe1_to_spe2_dir):
            os.makedirs(store_spe1_to_spe2_dir)

        ipt_spe2_syntenic_region_fl = input_output_dir + '/step01_build_ACR_in_syntenic_region_dir/' + '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir' + \
                                      '/temp_' + eachspe2 + '_syntenic_region_sorted.txt'

        ipt_spe1_syntenic_region_ACR_fl = input_output_dir + '/step01_build_ACR_in_syntenic_region_dir' + '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir' + \
                                      '/opt_intersect_' + s2_target_spe1_str + '_syn_ACR.txt'

        ##we will build a syntenic region ID file that only contains the ID line
        store_spe2_syntenic_region_dic = {}
        with open (ipt_spe2_syntenic_region_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                store_spe2_syntenic_region_dic[col[3]] = 1

        store_spe1_syntenic_region_ACR_dic = {}
        with open (ipt_spe1_syntenic_region_ACR_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                store_spe1_syntenic_region_ACR_dic[col[3]] = 1

        store_shared_syn_region_dic = {}
        for eachregion in store_spe2_syntenic_region_dic:
            if eachregion in store_spe1_syntenic_region_ACR_dic:
                store_shared_syn_region_dic[eachregion] = 1

        with open (store_spe1_to_spe2_dir + '/opt_shared_syntenic_region.txt','w+') as opt:
            for eachregion in store_shared_syn_region_dic:
                opt.write(eachregion + '\n')

        ipt_all_region_fl = store_spe1_to_spe2_dir + '/opt_shared_syntenic_region.txt'
        ipt_spe1_syn_region_ACR_fl = ipt_spe1_syntenic_region_ACR_fl
        ipt_spe2_syn_region_fl = ipt_spe2_syntenic_region_fl
        ipt_spe1_genome_fa_fl = input_all_spe_genome_dir + '/' + s2_target_spe1_str + '.fa'
        ipt_spe2_genome_fa_fl = input_all_spe_genome_dir + '/' + eachspe2 + '.fa'
        s2_functions.subfunction_run_parallel_blast (store_spe1_to_spe2_dir,ipt_all_region_fl,ipt_spe1_syn_region_ACR_fl,ipt_spe2_syn_region_fl,
                                    ipt_spe1_genome_fa_fl,ipt_spe2_genome_fa_fl,
                                    input_core_num)


def step03_merge_flt_res (input_output_dir,input_ct_acr_dir,s2_target_spe1_str,s2_target_spe2_str):

    step03_merge_flt_res_dir = input_output_dir + '/step03_merge_flt_res_dir'
    if not os.path.exists(step03_merge_flt_res_dir):
        os.makedirs(step03_merge_flt_res_dir)

    s2_target_spe2_str_list = s2_target_spe2_str.split(',')

    for eachspe2 in s2_target_spe2_str_list:

        store_spe1_to_spe2_dir = step03_merge_flt_res_dir + '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir'
        if not os.path.exists(store_spe1_to_spe2_dir):
            os.makedirs(store_spe1_to_spe2_dir)

        all_temp_opt_dir_list = glob.glob(input_output_dir + '/step02_conduct_blast_between_two_spe_dir/' + '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir' + '/temp_all_output_dir/*')

        all_blasted_fl_list = []
        for eachdir in all_temp_opt_dir_list:

            final_blasted_res_fl = eachdir + '/opt_final_blast_res.txt'
            all_blasted_fl_list.append(final_blasted_res_fl)

        cmd = 'cat ' + ' '.join(all_blasted_fl_list) + ' > ' + store_spe1_to_spe2_dir + '/opt_final_merged_blasted_res.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##filter the blasted res
        ##we will first make a new blast one with the same format as the FILTER SCRIPT

        ##we will first design the ACR ID
        ipt_ori_blast_fl = store_spe1_to_spe2_dir + '/opt_final_merged_blasted_res.txt'
        #ipt_spe1_ACR_fl = input_ACR_dir + '/' + s2_target_spe1_str + '_ACR.txt'
        ipt_spe2_syn_fl = input_output_dir + '/step01_build_ACR_in_syntenic_region_dir/' + \
                          '/store_' + s2_target_spe1_str + '_to_' + eachspe2 + '_dir/temp_' + eachspe2 + '_syntenic_region_sorted.txt'

        ##updating 053142
        if s2_target_spe1_str == 'Os':
            new_ct_spe1_str = 'Os'
        else:
            new_ct_spe1_str = s2_target_spe1_str

        ipt_spe1_ct_ACR_fl = input_ct_acr_dir + '/' + new_ct_spe1_str + '_acr_celltype.txt'


        s3_functions.modify_blasted_file (ipt_ori_blast_fl,ipt_spe1_ct_ACR_fl, ipt_spe2_syn_fl,store_spe1_to_spe2_dir)

        blast_file = store_spe1_to_spe2_dir + '/opt_final_blasted_modif_format.txt'
        evalue_threshold = 0.001
        output_file = store_spe1_to_spe2_dir + '/' + s2_target_spe1_str + '_ACRs.vs.' + eachspe2
        s3_functions.process_blast_file(blast_file, evalue_threshold, output_file)


def pipeline_analysis (input_ct_acr_dir,input_syntenic_region_fl,input_all_spe_genome_dir,input_output_dir,
                       input_core_num,open_step01,open_step02,input_spe1,input_spe2):

    target_spe1_str = input_spe1
    target_spe2_str_list = [input_spe2]

    if open_step01 == 'yes':

        store_all_spe_pair_res_dir = input_output_dir + '/store_all_spe_pair_res_dir'
        if not os.path.exists(store_all_spe_pair_res_dir):
            os.makedirs(store_all_spe_pair_res_dir)

        for eachspe2 in target_spe2_str_list:

            target_spe2_str = eachspe2

            store_opt_dir = store_all_spe_pair_res_dir + '/store_' + target_spe1_str + '_to_' + target_spe2_str + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            step01_build_ACR_in_syntenic_region(input_ct_acr_dir, input_syntenic_region_fl,
                                                target_spe1_str,
                                                target_spe2_str, store_opt_dir)

            step02_conduct_blast_between_two_spe(store_opt_dir, input_all_spe_genome_dir, target_spe1_str,
                                                 target_spe2_str, input_core_num)

            step03_merge_flt_res(store_opt_dir, input_ct_acr_dir, target_spe1_str, target_spe2_str)


    if open_step02 == 'yes':

        store_results_for_next_step_dir = input_output_dir + '/store_results_for_next_step_dir'
        if not os.path.exists(store_results_for_next_step_dir):
            os.makedirs(store_results_for_next_step_dir)



        collect_all_syntenic_regions_acrs_dir = store_results_for_next_step_dir + '/collect_all_syntenic_regions_acrs_dir'
        if not os.path.exists(collect_all_syntenic_regions_acrs_dir):
            os.makedirs(collect_all_syntenic_regions_acrs_dir)

        collect_output_dir_Pg_to_allspe_blast_dir = store_results_for_next_step_dir + '/collect_output_dir_Pg_to_allspe_blast_dir'
        if not os.path.exists(collect_output_dir_Pg_to_allspe_blast_dir):
            os.makedirs(collect_output_dir_Pg_to_allspe_blast_dir)

        collect_output_dir_Pg_to_allspe_blast_record_PgID_dir = store_results_for_next_step_dir + '/collect_output_dir_Pg_to_allspe_blast_record_PgID_dir'
        if not os.path.exists(collect_output_dir_Pg_to_allspe_blast_record_PgID_dir):
            os.makedirs(collect_output_dir_Pg_to_allspe_blast_record_PgID_dir)


        collect_all_other_spe_syntenic_genes_all_os_ACRs_dir = store_results_for_next_step_dir + '/collect_all_other_spe_syntenic_genes_all_os_ACRs_dir'
        if not os.path.exists(collect_all_other_spe_syntenic_genes_all_os_ACRs_dir):
            os.makedirs(collect_all_other_spe_syntenic_genes_all_os_ACRs_dir)


        for eachspe2 in target_spe2_str_list:

            spe1 = target_spe1_str
            spe2 = eachspe2

            ipt_spe1_to_spe2_dir = input_output_dir + '/store_all_spe_pair_res_dir/' + 'store_' + spe1 + '_to_' + spe2 + '_dir/'

            #################
            ##1) we will copy the Pg.syntenic_genes.all_Pg_ACRs.bed to the file 1
            opt_spe1_syntenic_gene_fl = ipt_spe1_to_spe2_dir  + '/step01_build_ACR_in_syntenic_region_dir/store_' + spe1 + '_to_' + spe2 + '_dir/' + spe1 + '.syntenic_genes.all_' + spe1 + '_ACRs.bed'

            cmd = 'cp ' + opt_spe1_syntenic_gene_fl + ' ' + store_results_for_next_step_dir + '/ipt_' + spe1 + '.syntenic_genes.ACR.bed'
            print(cmd)
            subprocess.call(cmd,shell=True)


            ##updating 073025
            ##we will store the ACR syntenic region first
            store_spe2_target_region_contain_ACR_dic = {}
            with open (store_results_for_next_step_dir + '/ipt_' + spe1 + '.syntenic_genes.ACR.bed','r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')
                    mt = re.match('(.+)_syntenic_.+',col[4])
                    spe2nm = mt.group(1)
                    if spe2nm == spe2:
                        target_region_ID = col[-1]
                        store_spe2_target_region_contain_ACR_dic[target_region_ID] = 1

            #################
            ##2) we will copy the syntenic region files to dir 2
            opt_syntenic_regions_spe2_fl = ipt_spe1_to_spe2_dir + '/step01_build_ACR_in_syntenic_region_dir/' + 'store_' + spe1 + '_to_' + spe2 + '_dir/' + 'temp_' + spe2 + '_syntenic_region_sorted.txt'
            ##/scratch/hy17471/database_031724/section3_add_pearlmillet_101424/10_downstream_analysis_071725/05_cross_species_072125/01_pipeline_blasted_ACRs_072125/output_dir/store_Pg_to_Uf_dir/step01_build_ACR_in_syntenic_region_dir/store_Pg_to_Uf_dir/temp_Uf_syntenic_region_sorted.txt
            store_final_line_list = []
            with open (opt_syntenic_regions_spe2_fl,'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')
                    regionID = col[3]
                    if regionID in store_spe2_target_region_contain_ACR_dic:
                        store_final_line_list.append(eachline)
            with open (collect_all_syntenic_regions_acrs_dir + '/' +  spe2 + '.syntenic_regions.os_acrs.bed','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            #cmd = 'cp ' + opt_syntenic_regions_spe2_fl + ' ' + collect_all_syntenic_regions_acrs_dir + '/' + spe2 + '.syntenic_regions.os_acrs.bed'
            #print(cmd)
            #subprocess.call(cmd,shell=True)

            opt_syntenic_regions_spe1_fl = ipt_spe1_to_spe2_dir + '/step01_build_ACR_in_syntenic_region_dir/' + 'store_' + spe1 + '_to_' + spe2 + '_dir/' + 'temp_' + spe1 + '_syntenic_region_sorted.txt'
            cmd = 'cp ' + opt_syntenic_regions_spe1_fl + ' ' + collect_all_syntenic_regions_acrs_dir + '/' + spe1 + '.syntenic_regions.os_acrs.bed'
            print(cmd)
            subprocess.call(cmd, shell=True)

            ################
            ##3) we will copy the all blast dir
            opt_blasted_fl = ipt_spe1_to_spe2_dir + '/step03_merge_flt_res_dir/' + 'store_' + spe1 + '_to_' + spe2 + '_dir/' + spe1 + '_ACRs.vs.' + spe2 + '.blast_passing_regions.intersecting_regions.bed'
            cmd = 'cp ' + opt_blasted_fl + ' ' + collect_output_dir_Pg_to_allspe_blast_dir + '/' + spe1 + '_ACRs.vs.' + spe2 + '_blast.bed'
            print(cmd)
            subprocess.call(cmd,shell=True)

            ################
            ##4) we will copy the all blast ref file
            opt_blasted_ref_fl =  ipt_spe1_to_spe2_dir + '/step03_merge_flt_res_dir/' + 'store_' + spe1 + '_to_' + spe2 + '_dir/' + spe1 + '_ACRs.vs.' + spe2 + '.blast_passing_regions.intersecting_regions.ref.bed'
            cmd = 'cp ' + opt_blasted_ref_fl + ' ' + collect_output_dir_Pg_to_allspe_blast_record_PgID_dir + '/' + spe1 + '_ACRs.vs.' + spe2 + '_blast_record_' + spe1 + 'ID.bed'
            print(cmd)
            subprocess.call(cmd,shell=True)


            #################
            ##5) we will build the spe1 syntenic gene file
            ##it is in the input_syntenic_region_fl
            ##we will generate the gene location within each location in this file
            ##and we will copy them into the collect_all_other_spe_syntenic_genes_all_os_ACRs_dir
            store_target_spe_line_list = []
            count = 0
            with open (input_syntenic_region_fl,'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split('\t')
                    count += 1
                    if count != 1:
                        spe2_chr = col[4]
                        spe2_st = col[5]
                        spe2_ed = col[6]
                        spe2_geneID = col[7]
                        region_ID = col[-2]
                        spenm = col[8]

                        if spenm == spe2:

                            #updating 073025
                            if region_ID in store_spe2_target_region_contain_ACR_dic:

                                final_line = spe2_chr + '\t' + spe2_st + '\t' + spe2_ed + '\t' + spe2_geneID + '\t' + region_ID
                                store_target_spe_line_list.append(final_line)

            with open (collect_all_other_spe_syntenic_genes_all_os_ACRs_dir + '/' + spe2 + '.syntenic_genes.all_' + spe2 + '_ACRs.bed','w+') as opt:
                for eachline in store_target_spe_line_list:
                    opt.write(eachline + '\n')





pipeline_analysis (input_ct_acr_dir,input_syntenic_region_fl,input_all_spe_genome_dir,input_output_dir,
                       input_core_num,open_step01,open_step02,input_spe1,input_spe2)













