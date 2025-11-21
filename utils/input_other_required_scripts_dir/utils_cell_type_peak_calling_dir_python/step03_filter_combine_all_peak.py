#!/usr/bin/env python

import re
import glob
import subprocess
import os

##updating 112125 we will remove the previous bed file
##updating 062425 add a function to filter out peaks that is out of genome size per chromosome



def filter_acr (input_output_dir,input_fdr_Rscript,finalfdr,new_fdr_dir_output_dir):

    store_control_peak_intersect_dir = input_output_dir + '/store_control_peak_intersect_dir'
    store_raw_peak_intersect_dir = input_output_dir + '/store_raw_peak_intersect_dir'

    store_filtered_peak_dir = new_fdr_dir_output_dir + '/store_filtered_peak_dir'
    if not os.path.exists(store_filtered_peak_dir):
        os.makedirs(store_filtered_peak_dir)

    raw_peak_fl_list = glob.glob(store_raw_peak_intersect_dir + '/*')

    for eachfl in raw_peak_fl_list:

        mt = re.match('.+/(.+)',eachfl)
        flnm = mt.group(1)

        mt = re.match('cluster\.(.+)\.macs2_peaks\.flt\.narrowPeak',flnm)

        organ_ct = mt.group(1)

        control_peak_fl = store_control_peak_intersect_dir + '/opt_' + organ_ct + '_permute_control_peak_sorted_intersect_tn5.txt'

        cmd = 'Rscript ' + input_fdr_Rscript + \
              ' ' + eachfl + \
              ' ' + control_peak_fl + \
              ' ' + finalfdr + \
              ' ' + organ_ct + \
              ' ' + store_filtered_peak_dir
        subprocess.call(cmd,shell=True)

def check_number_each_celltype (input_output_dir,input_meta_fl,new_fdr_dir_output_dir,
                                s2_targettn5_colNum,target_colnm):

    ##updating 031625
    ##some times the annotation is not the final
    store_organ_ct_tn5_dic = {}
    store_cellnum_dic = {}
    count = 0
    target_index = 0
    with open (input_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1

            if count == 1:
                temp_target_index = col.index(target_colnm)

                all_len_col = len(col)

                target_index = temp_target_index - all_len_col

            else:

                ##updating 111622
                #organ_ct = col[int(s1_targetclust_colNum) - 1]
                organ_ct = col[int(target_index)]
                #organ_ct = col[24]
                tn5 = col[int(s2_targettn5_colNum)]

                if organ_ct in store_cellnum_dic:
                    store_cellnum_dic[organ_ct] += 1
                else:
                    store_cellnum_dic[organ_ct] = 1

                if organ_ct in store_organ_ct_tn5_dic:
                    store_organ_ct_tn5_dic[organ_ct] += int(tn5)
                else:
                    store_organ_ct_tn5_dic[organ_ct] = int(tn5)


    store_filtered_peak_dir = new_fdr_dir_output_dir + '/store_filtered_peak_dir'

    store_raw_peak_dir = input_output_dir + '/store_MACS2_raw_peak_dir'

    store_final_line_list = []
    peak_fl_list = glob.glob(store_filtered_peak_dir + '/*')
    for eachpeakfl in peak_fl_list:

        mt = re.match('.+/(.+)',eachpeakfl)
        flnm = mt.group(1)

        mt = re.match('opt_flt_(.+)_peak.txt',flnm)
        organ_ct = mt.group(1)

        flt_number_peak = 0
        with open (eachpeakfl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                flt_number_peak += 1

        raw_peak_fl = store_raw_peak_dir + '/cluster.' + organ_ct + '.macs2_peaks.flt.narrowPeak'
        raw_number_peak = 0
        with open (raw_peak_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                raw_number_peak += 1

        final_line = organ_ct + '\t' + str(store_organ_ct_tn5_dic[organ_ct]) + '\t' + str(store_cellnum_dic[organ_ct]) + '\t' + str(flt_number_peak) + '\t' + str(raw_number_peak)
        store_final_line_list.append(final_line)

    with open (new_fdr_dir_output_dir + '/opt_summary_peak_num_for_eachCT.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


##uupdating 040122 combine the peak
def combine_peak_from_clusters (input_output_dir,input_required_script_dir,prefix,new_fdr_dir_output_dir):

    store_summits_peak_dir = new_fdr_dir_output_dir + '/store_summits_peak_dir'
    if not os.path.exists(store_summits_peak_dir):
        os.makedirs(store_summits_peak_dir)

    store_filtered_summits_peak_dir = new_fdr_dir_output_dir + '/store_filtered_summits_peak_dir'
    if not os.path.exists(store_filtered_summits_peak_dir):
        os.makedirs(store_filtered_summits_peak_dir)

    store_final_unique_peak_dir = new_fdr_dir_output_dir + '/store_final_unique_peak_dir'
    if not os.path.exists(store_final_unique_peak_dir):
        os.makedirs(store_final_unique_peak_dir)

    store_MACS2_calling_dir = input_output_dir + '/store_MACS2_calling_dir'

    store_filtered_peak_dir = new_fdr_dir_output_dir + '/store_filtered_peak_dir'

    all_temp_dir_list = glob.glob(store_MACS2_calling_dir + '/*')

    for eachdir in all_temp_dir_list:

        cmd = 'cp ' + eachdir + '/*_summits.bed ' + store_summits_peak_dir
        subprocess.call(cmd,shell=True)

    ##now get summit of the filtered peaks
    all_flt_fl_list = glob.glob(store_filtered_peak_dir + '/*')

    for eachfl in all_flt_fl_list:

        mt = re.match('.+/(.+)\.txt',eachfl)
        flnm = mt.group(1)

        mt = re.match('opt_flt_(.+)_peak',flnm)
        organ_celltype = mt.group(1)

        opt_submmit_peak_fl = store_summits_peak_dir + '/cluster.' + organ_celltype + '.macs2_summits.bed'

        ##check the summit
        cmd = 'bedtools intersect -a ' + opt_submmit_peak_fl + ' -b ' + eachfl + ' -u > ' + \
              store_filtered_summits_peak_dir + '/cluster.' + organ_celltype + '.flt_summits.bed'
        subprocess.call(cmd,shell=True)

    ##allow the script to be the current working dir
    ##adjust the peak
    ori_cwd = os.getcwd()

    adjust_peak_script = input_required_script_dir + '/adjustPeaks.sh'
    normalize_score_script = input_required_script_dir + '/normalize_score.py'
    selectNonOverlapping_script = input_required_script_dir + '/selectNonOverlapping.py'

    ##udpating 061825
    expandannotatepeak_script = input_required_script_dir + '/expand_and_annotate_peaks.py'

    cmd = 'cp ' + adjust_peak_script + ' ' + normalize_score_script + ' ' + selectNonOverlapping_script + ' ' + expandannotatepeak_script + ' ' + store_filtered_summits_peak_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##open the store_peak_calling_dir
    os.chdir(store_filtered_summits_peak_dir)
    #cmd = 'cd ' + store_peak_calling_dir
    #subprocess.call(cmd,shell=True)

    ##now current diretory is store_peak_calling_dir
    ##call peaks
    cmd = 'bash adjustPeaks.sh' + \
          ' ' + prefix
    subprocess.call(cmd,shell=True)

    ##go back to the originial cwd
    #cmd = 'cd ' + ori_cwd
    #subprocess.call(cmd,shell=True)
    os.chdir(ori_cwd)

    cmd = 'cp ' + store_filtered_summits_peak_dir + '/' + prefix + '.unique500bpPeaks.bed ' + store_final_unique_peak_dir
    subprocess.call(cmd,shell=True)



def FDR_filteration (input_step02_output_dir,utils_cell_type_peak_calling_dir,
                     input_meta_fl,finalfdr_list,s2_targettn5_colNum,s3_open_filter_peak_final_dir,pval_or_qval,organ_pqval_num,target_colnm):

    for eachFDRnum in finalfdr_list:

        target_organ_nm_pqval_FDRnum_dir = input_step02_output_dir + '/' + eachFDRnum
        if not os.path.exists(target_organ_nm_pqval_FDRnum_dir):
            os.makedirs(target_organ_nm_pqval_FDRnum_dir)

        input_fdr_Rscript = utils_cell_type_peak_calling_dir + '/filter_ACRs_eFDR.R'

        ##we would first load the target_organ_nm_pqval_dir since some inputs would be in this directory
        filter_acr(input_step02_output_dir, input_fdr_Rscript, eachFDRnum, target_organ_nm_pqval_FDRnum_dir)
        check_number_each_celltype(input_step02_output_dir, input_meta_fl, target_organ_nm_pqval_FDRnum_dir,
                                   s2_targettn5_colNum,target_colnm)
        combine_peak_from_clusters(input_step02_output_dir, utils_cell_type_peak_calling_dir, 'final',
                                   target_organ_nm_pqval_FDRnum_dir)

        ##collect all the peak information
        cmd = 'cp ' + target_organ_nm_pqval_FDRnum_dir + '/store_final_unique_peak_dir/' + 'final' + '.unique500bpPeaks.bed ' + \
              s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks.bed'
        print(cmd)
        subprocess.call(cmd, shell=True)


        ##updating  122824
        ##change the last col to be ACR number
        store_final_line_list = []
        count = 0
        with open (s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks.bed','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')
                count += 1
                ACRnm = 'ACR_' + str(count)
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + ACRnm
                store_final_line_list.append(final_line)

        with open (s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks_new.bed','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        cmd = 'rm ' + s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks.bed'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'mv ' + s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks_new.bed ' + s3_open_filter_peak_final_dir + '/opt_' + 'final' + '.' + pval_or_qval + organ_pqval_num + '.FDR' + eachFDRnum + '.unique500bpPeaks.bed'
        print(cmd)
        subprocess.call(cmd,shell=True)




##updating 062425
def filter_peak_by_chr_size (ipt_filtered_peak_fl, ipt_fa_fai_fl,ipt_opt_dir):


    mt = re.match('.+/(.+)\.bed',ipt_filtered_peak_fl)
    flnm = mt.group(1)


    store_chr_max_size_dic = {}
    with open (ipt_fa_fai_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            chrnm = col[0]
            chrsize = col[1]
            store_chr_max_size_dic[chrnm] = chrsize

    store_chr_ACRloclist_dic = {}
    with open (ipt_filtered_peak_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            chrnm = col[0]
            acrloc = col[0] + '\t' + col[1] + '\t' + col[2]

            if chrnm in store_chr_ACRloclist_dic:
                store_chr_ACRloclist_dic[chrnm].append(acrloc)
            else:
                store_chr_ACRloclist_dic[chrnm] = []
                store_chr_ACRloclist_dic[chrnm].append(acrloc)

    store_final_line_list = []
    acrid = 0
    for eachchrnm in store_chr_max_size_dic:
        chrmax_size = store_chr_max_size_dic[eachchrnm]

        if eachchrnm in store_chr_ACRloclist_dic:

            acrloclist = store_chr_ACRloclist_dic[eachchrnm]

            count = 0
            for eachacrloc in acrloclist:

                acrid += 1

                count += 1
                if count == len(acrloclist):

                    lastacr = eachacrloc
                    lastacr_list = lastacr.split('\t')
                    end_loc = lastacr_list[-1]

                    if int(end_loc) < int(chrmax_size):

                        ##we will keep this line
                        final_acr = lastacr + '\t' + str('ACR_' + str(acrid))
                        store_final_line_list.append(final_acr)

                else:
                    final_acr = eachacrloc + '\t' + str('ACR_' + str(acrid))
                    store_final_line_list.append(final_acr)


    with open (ipt_opt_dir  + '/temp_' + flnm + '.fixsize.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + ipt_opt_dir  + '/temp_' + flnm + '.fixsize.bed > ' + \
          ipt_opt_dir + '/' + flnm + '.fx.sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'rm ' + ipt_opt_dir  + '/temp_' + flnm + '.fixsize.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)


    ##updating 112125 remove the previous file and only retain the
    cmd = 'rm ' + ipt_opt_dir + '/' + flnm + '.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)










