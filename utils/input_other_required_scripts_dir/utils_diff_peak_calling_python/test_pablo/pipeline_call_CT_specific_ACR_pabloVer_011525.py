#!/usr/bin/env python


import re
import glob
import sys
import subprocess
import os

input_peak_tn5_fl = sys.argv[1]
##/public2/home/yanhaidong/working_dir/yanhaidong/Socrates2/call_diff_peaks/output_dir/s1_open_prepare_peak_tn5_final/opt_peak_sparse_sorted_dir/opt_peak_sorted.sparse

##/public2/home/yanhaidong/working_dir/yanhaidong/Socrates2/call_diff_peaks/output_dir/s1_open_prepare_peak_tn5_final/opt_peak_sparse_sorted_dir/opt_peak_sorted.sparse

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_5_iden_DApeaks_040422/01_make_peak_sparse_fl_040422/output_dir_051323/opt_all_peak_sparse_nineorgans_sorted.sparse

input_meta_fl = sys.argv[2]
##/public2/home/yanhaidong/working_dir/yanhaidong/Socrates2/marker_annot/marker_scAClass_output/opt_meta_ScAClass_predicted_cell_identity.txt

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/opt_final_comb_annot_ver8.4.txt

input_final_peak_fl = sys.argv[3]
##/public2/home/yanhaidong/working_dir/yanhaidong/Socrates2/call_peaks_per_cell_type/output_dir/s3_open_filter_peak_final_dir/opt_final.pval0.5.FDR0.005.unique500bpPeaks_choose.bed

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_6_characterize_ACRs_060322/output_dir_v8.2_update_051323/step01_ACR_summary/s1_loc_analysis/opt_ACR_sorted_rmdup.txt

input_ctsACR_R_script_fl = sys.argv[4]
#


input_configure_fl = sys.argv[5]


input_output_dir = sys.argv[6]



def step01_prepare_peak_acc (input_peak_tn5_fl,input_final_peak_fl,input_output_dir):

    step01_prepare_peak_acc_dir = input_output_dir + '/step01_prepare_peak_acc_dir'
    if not os.path.exists(step01_prepare_peak_acc_dir):
        os.makedirs(step01_prepare_peak_acc_dir)

    store_final_peak_dic = {}
    with open (input_final_peak_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            final_line = col[0] + '_' + col[1] + '_' + col[2]
            store_final_peak_dic[final_line] = 1

    store_all_acr_dic = {}
    with open (input_peak_tn5_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            new_peak = col[0].replace('chr','')
            acrloc = '\t'.join(new_peak.split('_'))
            store_all_acr_dic[acrloc] = new_peak

    store_acr_ID_dic = {}
    store_final_line_list = []
    count = 0
    for eachacrloc in store_all_acr_dic:

        acrloc_nospace = store_all_acr_dic[eachacrloc]

        count += 1
        acrID = 'scACR_' + str(count)

        if acrloc_nospace in store_final_peak_dic:

            final_line = eachacrloc + '\t' + acrID + '\t' + '1'
            store_final_line_list.append(final_line)

            store_acr_ID_dic[acrloc_nospace] = acrID

    with open (step01_prepare_peak_acc_dir + '/opt_peak.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    tn5_fl = open (input_peak_tn5_fl,'r')
    opt_fl = open (step01_prepare_peak_acc_dir + '/opt_accessibility.txt','w')

    for eachline in tn5_fl:

        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        acrlocnospace = col[0]
        acrlocnospace = acrlocnospace.replace('chr','')

        if acrlocnospace in store_final_peak_dic:

            acrID = store_acr_ID_dic[acrlocnospace]

            final_line = acrID + '\t' + col[1] + '\t' + col[2]
            opt_fl.write(final_line + '\n')

    tn5_fl.close()
    opt_fl.close()

def step02_call_CT_PabloVer (input_ctsACR_R_script_fl,input_meta_fl,input_output_dir,
                             s2_prefix,s2_target_cluster,
                             s2_threshold,s2_stat_test,s2_null_permutations,s2_entropy_bootstraps):

    step02_call_CT_PabloVer_dir = input_output_dir + '/step02_call_CT_PabloVer_dir'
    if not os.path.exists(step02_call_CT_PabloVer_dir):
        os.makedirs(step02_call_CT_PabloVer_dir)

    ipt_acr_acc_fl = input_output_dir + '/step01_prepare_peak_acc_dir/opt_accessibility.txt'
    ipt_acr_peak_fl = input_output_dir + '/step01_prepare_peak_acc_dir/opt_peak.bed'

    cmd = 'Rscript ' + input_ctsACR_R_script_fl + \
          ' --input_data ' + ipt_acr_acc_fl + \
          ' --peak_file ' + ipt_acr_peak_fl + \
          ' --meta ' + input_meta_fl + \
          ' --meta_slot ' + s2_target_cluster+ \
          ' --prefix ' + s2_prefix + \
          ' --stat_test ' + s2_stat_test + \
          ' --null_permutations ' + s2_null_permutations + \
          ' --entropy_bootstraps ' + s2_entropy_bootstraps
    print(cmd)
    subprocess.call(cmd,shell=True)




store_target_parameter_dic = {}
with open (input_configure_fl,'r') as ipt:
    for eachline in ipt:
        eachline = eachline.strip('\n')
        if not eachline.startswith('#'):
            col = eachline.strip().split('=')
            store_target_parameter_dic[col[0]] = col[1]

step01 = store_target_parameter_dic['step01'] ##yes or no
step02 = store_target_parameter_dic['step02'] ##yes or no

s2_prefix = store_target_parameter_dic['s2_prefix']
s2_target_cluster = store_target_parameter_dic['s2_target_cluster']
s2_threshold = store_target_parameter_dic['s2_threshold']
s2_stat_test = store_target_parameter_dic['s2_stat_test']
s2_null_permutations = store_target_parameter_dic['s2_null_permutations']
s2_entropy_bootstraps = store_target_parameter_dic['s2_entropy_bootstraps']

if step01 == 'yes':

    step01_prepare_peak_acc (input_peak_tn5_fl,input_final_peak_fl,input_output_dir)

if step02 == 'yes':
    step02_call_CT_PabloVer(input_ctsACR_R_script_fl, input_meta_fl, input_output_dir,
                            s2_prefix, s2_target_cluster,
                            s2_threshold, s2_stat_test, s2_null_permutations, s2_entropy_bootstraps)



















