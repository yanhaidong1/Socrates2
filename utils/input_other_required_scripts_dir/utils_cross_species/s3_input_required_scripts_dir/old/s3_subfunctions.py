#!/usr/bin/env python



import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats


def subfunction_check_all_maize_region_InOrNotH3K27me3 (ipt_target_organ_H3K27_peak_spe1_fl,ipt_target_organ_H3K27_peak_spe2_fl,
                                                        ipt_final_summary_blast_fl,
                                                        ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                        opt_dir,
                                                        s3_extend_H3K27me3_flank_bp,ipt_spe1_prefix,ipt_spe2_prefix):

    ##updating 111023
    ##add the rice gene cate
    store_acr_genecate_dist_spe1_dic = {}
    with open (ipt_spe1_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_acr_genecate_dist_spe1_dic[col[0]] = {'cate':col[1],'dist':col[2]}

    store_acr_genecate_dist_spe2_dic = {}
    with open (ipt_spe2_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_acr_genecate_dist_spe2_dic[col[0]] = {'cate':col[1],'dist':col[2]}


    store_maize_loc_dic = {}
    store_maize_acr_loc_dic = {}
    store_rice_acr_loc_dic = {}
    count = 0
    with open (ipt_final_summary_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                maize_region = col[4]
                if maize_region != 'none':
                    maize_loc_line = '\t'.join(maize_region.split('_'))
                    store_maize_loc_dic[maize_loc_line] = 1

                    maize_acr = col[5]
                    if maize_acr != 'none':
                        maize_acr_loc_line = '\t'.join(maize_acr.split('_'))
                        store_maize_acr_loc_dic[maize_acr_loc_line] = 1

                rice_acr = col[1]
                if rice_acr != 'none':
                    rice_acr_loc_line = '\t'.join(rice_acr.split('_'))
                    store_rice_acr_loc_dic[rice_acr_loc_line] = 1



    with open (opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region.txt','w+') as opt:
        for eachloc in store_maize_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region.txt > ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    with open (opt_dir + '/temp_' + ipt_spe2_prefix + '_acr.txt','w+') as opt:
        for eachloc in store_maize_acr_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr.txt' + ' > ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    with open (opt_dir + '/temp_' + ipt_spe1_prefix + '_acr.txt','w+') as opt:
        for eachloc in store_rice_acr_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr.txt > ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    with open(ipt_target_organ_H3K27_peak_spe2_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            st = int(col[1]) - int(s3_extend_H3K27me3_flank_bp)
            if st < 0:
                new_st = 0
            else:
                new_st = st
            ed = int(col[2]) + int(s3_extend_H3K27me3_flank_bp)

            final_line = col[0] + '\t' + str(new_st) + '\t' + str(ed)
            store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak.txt',
            'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak.txt' + ' > ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##intersection
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt > ' + \
          opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_loc_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_maize_loc_blastedToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_loc_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            blastregion_loc = col[0] + '_' + col[1] + '_' + col[2]
            store_maize_loc_blastedToH3K27me3_dic[blastregion_loc] = 1

    ##intersect with acr
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr_sorted.txt' + \
          ' -b ' +  opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt > ' + \
           opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_acr_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_maize_acrToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_acr_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            maize_acr = col[0] + '_' + col[1] + '_' + col[2]
            store_maize_acrToH3K27me3_dic[maize_acr] = 1


    ##intersect with acr in the rice
    store_final_line_list = []
    with open (ipt_target_organ_H3K27_peak_spe1_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            st = int(col[1]) - int(s3_extend_H3K27me3_flank_bp)
            if st < 0:
                new_st = 0
            else:
                new_st = st
            ed = int(col[2]) + int(s3_extend_H3K27me3_flank_bp)

            final_line = col[0] + '\t' + str(new_st) + '\t' + str(ed)
            store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_' + ipt_spe1_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_' + ipt_spe1_prefix + '.txt > ' + \
          opt_dir + '/temp_intersect_' + ipt_spe1_prefix + '_acr_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_rice_acrToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe1_prefix + '_acr_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            rice_acr = col[0] + '_' + col[1] + '_' + col[2]
            store_rice_acrToH3K27me3_dic[rice_acr] = 1


    ##updating 110723
    ##Here we will build a new file updating if the ACR overlapping with the H3K27me3
    store_final_line_list = []
    count = 0
    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACR = col[1]
                maizeACR = col[5]
                maizeregion = col[4]
                if riceACR in store_rice_acrToH3K27me3_dic:
                    riceOverH3K27Cate = 'OverlapH3K27me3'
                else:
                    riceOverH3K27Cate = 'none'

                if maizeACR in store_maize_acrToH3K27me3_dic:
                    maizeOVerH3K27Cate = 'OverlapH3K27me3'
                else:
                    maizeOVerH3K27Cate = 'none'

                if riceACR in store_acr_genecate_dist_spe1_dic:
                    rice_genecate = store_acr_genecate_dist_spe1_dic[riceACR]['cate']
                    rice_genedist = store_acr_genecate_dist_spe1_dic[riceACR]['dist']
                else:
                    rice_genecate = 'none'
                    rice_genedist = 'none'

                if maizeACR in store_acr_genecate_dist_spe2_dic:
                    maize_genecate = store_acr_genecate_dist_spe2_dic[maizeACR]['cate']
                    maize_genedist = store_acr_genecate_dist_spe2_dic[maizeACR]['dist']

                else:
                    maize_genecate = 'none'
                    maize_genedist = 'none'

                if maizeregion in store_maize_loc_blastedToH3K27me3_dic:
                    maizeregionOverH3K27Cate = 'OverlapH3K27me3'
                else:
                    maizeregionOverH3K27Cate = 'none'

                final_line = eachline + '\t' + riceOverH3K27Cate + '\t' + maizeOVerH3K27Cate + '\t' + rice_genecate + '\t' + maize_genecate + '\t' + rice_genedist + '\t' + maize_genedist + '\t' + maizeregionOverH3K27Cate
                store_final_line_list.append(final_line)
            else:
                final_line = eachline + '\t' + ipt_spe1_prefix + 'OverH3K27Cate' + '\t' + ipt_spe2_prefix + 'OVerH3K27Cate' + '\t' + ipt_spe1_prefix + 'GeneCate' + '\t' + ipt_spe2_prefix + 'GeneCate' + \
                            '\t' + ipt_spe1_prefix + 'GeneDist' + '\t' + ipt_spe2_prefix + 'GeneDist' + '\t' + ipt_spe2_prefix + 'RegionOverH3K27Cate'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_summary_file_add_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    return (store_maize_loc_blastedToH3K27me3_dic,store_maize_acrToH3K27me3_dic,store_rice_acrToH3K27me3_dic)


def subfunction_store_different_cate_riceACRs (maizeACRloc,maize_H3K27me3Cate,maize_BroOrRes,store_maize_acrToH3K27me3_dic,
                                               riceACRloc,maizeblastedRg,store_maize_loc_blastedToH3K27me3_dic,
                                               store_broad_to_broad_riceACRnotCoverH3K27me3_dic,store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic,
                                               store_broad_to_noACR_notCoverH3K27me3_dic):



    ##There are something wrong for some

    ##for the broad to broad we will consider two cases under the H3K27me3 and not for the maize
    if maizeACRloc != 'none':

        if maize_H3K27me3Cate == 'BroadInflankH3K27me3peak':
            cate = 'H3K27me3'

            if cate in store_broad_to_broad_riceACRnotCoverH3K27me3_dic:
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
            else:
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate] = {}
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

        else:

            if maize_H3K27me3Cate != 'RestrictedInflankH3K27me3peak':

                ##make sure maizeACRloc not cover H3K27me3
                if maizeACRloc not in store_maize_acrToH3K27me3_dic:

                    if maize_BroOrRes == 'broad':
                        cate = 'H3K27me3Absent'

                        if cate in store_broad_to_broad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                    else:

                        cate = 'H3K27me3Absent'
                        if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                ##another cate to include the others underlying the H3K27me3
                else:

                    if maize_BroOrRes == 'others' or maize_BroOrRes == 'restricted':
                        cate = 'H3K27me3'
                        if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                    ##Here the else the broad we already consider the above condition


            else:
                ##Here is only for the restriced in H3K27me3peak we will also check the others in the H3K27me3
                cate = 'H3K27me3'
                if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1


    ## if there is no acr
    else:
        if maizeblastedRg != 'none':

            ##for the H3K27me3
            if maizeblastedRg in store_maize_loc_blastedToH3K27me3_dic:
                cate = 'H3K27me3'
                if cate in store_broad_to_noACR_notCoverH3K27me3_dic:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate] = {}
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
            else:
                cate = 'H3K27me3Absent'
                if cate in store_broad_to_noACR_notCoverH3K27me3_dic:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate] = {}
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1



def subfunction_summarize_overview_for_H3K27me3 (ipt_final_summary_blast_fl,ipt_target_organ_H3K27_peak_spe1_fl,
                                                 ipt_target_organ_H3K27_peak_spe2_fl,
                                                 ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                 spe1_prefix,spe2_prefix,opt_dir,
                                                 s3_extend_H3K27me3_flank_bp):

    ##we need to first store blasted loc whether they are overlapped with the H3K27me3
    store_maize_loc_blastedToH3K27me3_dic, store_maize_acrToH3K27me3_dic,\
    store_rice_acrToH3K27me3_dic = subfunction_check_all_maize_region_InOrNotH3K27me3(ipt_target_organ_H3K27_peak_spe1_fl,
                                                                                      ipt_target_organ_H3K27_peak_spe2_fl,
                                                                                      ipt_final_summary_blast_fl,
                                                                                      ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                                                      opt_dir,
                                                                                      s3_extend_H3K27me3_flank_bp, spe1_prefix,spe2_prefix)

    #################
    ##category type 1

    ##H3K27me3 not shared
    count = 0
    store_H3K27me3_absent_riceACR_dic = {}
    store_H3K27me3_sharedbroad_acc_riceACR_dic = {}
    store_H3K27me3_sharedrestricted_acc_riceACR_dic = {}
    store_H3K27me3_sharedinacc_riceACR_dic = {}
    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRloc = col[1]
                maizeACRloc = col[5]
                maizeblastedRg = col[4]
                rice_H3K27me3Cate = col[2]
                maize_H3K27me3Cate = col[6]

                if rice_H3K27me3Cate == 'BroadInflankH3K27me3peak':

                    if maizeACRloc != 'none':
                        ##indicate maize ACR not overlap with the H3K27me3
                        if maizeACRloc not in store_maize_acrToH3K27me3_dic:
                            store_H3K27me3_absent_riceACR_dic[riceACRloc] = 1

                        else:
                            if maize_H3K27me3Cate == 'BroadInflankH3K27me3peak':
                                store_H3K27me3_sharedbroad_acc_riceACR_dic[riceACRloc] = 1

                            if maize_H3K27me3Cate == 'RestrictedInflankH3K27me3peak':
                                store_H3K27me3_sharedrestricted_acc_riceACR_dic[riceACRloc] = 1

                    else:
                        if maizeblastedRg != 'none':
                            if maizeblastedRg not in store_maize_loc_blastedToH3K27me3_dic:
                                store_H3K27me3_absent_riceACR_dic[riceACRloc] = 1
                            else:
                                store_H3K27me3_sharedinacc_riceACR_dic[riceACRloc] = 1

    ##report the results
    store_final_line_list = []
    final_line = 'H3K27me3_absent' + '\t' + str(len(list(store_H3K27me3_absent_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_broad_accessible' + '\t' + str(len(list(store_H3K27me3_sharedbroad_acc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_restricted_accessible' + '\t' + str(len(list(store_H3K27me3_sharedrestricted_acc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_inaccessible' + '\t' + str(len(list(store_H3K27me3_sharedinacc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    #################
    ##category type 2
    ##This is for the covering H3K27me3
    store_broad_to_broad_riceACR_dic = {}
    store_broad_to_notbroad_riceACR_dic = {}
    store_broad_to_noACR_dic = {}

    ##updating 110623
    ##Here we will build a background to use the fisher to check the correspondings were enriched for the ACR under the H3K27me3
    store_broad_to_broad_riceACRnotCoverH3K27me3_dic = {}
    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic = {}
    store_broad_to_noACR_notCoverH3K27me3_dic = {}

    store_notbroad_to_broad_riceACR_dic = {}
    store_notbroad_to_notbroad_riceACR_dic = {}
    store_notbroad_to_noACR_dic = {}

    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRloc = col[1]
                rice_BroOrRes = col[3]
                maizeACRloc = col[5]
                maizeblastedRg = col[4]
                rice_H3K27me3Cate = col[2]
                maize_H3K27me3Cate = col[6]
                maize_BroOrRes = col[7]

                if rice_H3K27me3Cate == 'BroadInflankH3K27me3peak':

                    subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                              store_maize_acrToH3K27me3_dic,
                                                              riceACRloc, maizeblastedRg,
                                                              store_maize_loc_blastedToH3K27me3_dic,
                                                              store_broad_to_broad_riceACR_dic,
                                                              store_broad_to_notbroad_riceACR_dic,
                                                              store_broad_to_noACR_dic)


                #####################################################
                ##Here we will store the rice ACR not in the H3K27me3
                else:
                    if riceACRloc not in store_rice_acrToH3K27me3_dic:

                        if rice_BroOrRes == 'broad':

                            subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                                      store_maize_acrToH3K27me3_dic,
                                                                      riceACRloc, maizeblastedRg,
                                                                      store_maize_loc_blastedToH3K27me3_dic,
                                                                      store_broad_to_broad_riceACRnotCoverH3K27me3_dic,
                                                                      store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic,
                                                                      store_broad_to_noACR_notCoverH3K27me3_dic)


                    else:


                        subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                                  store_maize_acrToH3K27me3_dic,
                                                                  riceACRloc, maizeblastedRg,
                                                                  store_maize_loc_blastedToH3K27me3_dic,
                                                                  store_notbroad_to_broad_riceACR_dic,
                                                                  store_notbroad_to_notbroad_riceACR_dic,
                                                                  store_notbroad_to_noACR_dic)


    ##now we will check the proportion for all cate
    store_final_line_list = []

    total_broad_to_broad_riceACR_dic = {}
    for eachcate in store_broad_to_broad_riceACR_dic:
        for eachacr in store_broad_to_broad_riceACR_dic[eachcate]:
            total_broad_to_broad_riceACR_dic[eachacr] = 1
    total_broad_to_broad_riceACR_num = len(list(total_broad_to_broad_riceACR_dic.keys()))
    for eachcate in store_broad_to_broad_riceACR_dic:
        cate_acr_dic = store_broad_to_broad_riceACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToBroad' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_broad_riceACR_num)
        store_final_line_list.append(final_line)

    total_broad_to_notbroad_riceACR_dic = {}
    for eachcate in store_broad_to_notbroad_riceACR_dic:
        for eachacr in store_broad_to_notbroad_riceACR_dic[eachcate]:
            total_broad_to_notbroad_riceACR_dic[eachacr] = 1
    total_broad_to_notbroad_riceACR_num = len(list(total_broad_to_notbroad_riceACR_dic.keys()))
    for eachcate in store_broad_to_notbroad_riceACR_dic:
        cate_acr_dic = store_broad_to_notbroad_riceACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToNotBroad' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_notbroad_riceACR_num)
        store_final_line_list.append(final_line)

    total_broad_to_noACR_riceACR_dic = {}
    for eachcate in store_broad_to_noACR_dic:
        for eachacr in store_broad_to_noACR_dic[eachcate]:
            total_broad_to_noACR_riceACR_dic[eachacr] = 1
    total_broad_to_noACR_riceACR_num = len(list(total_broad_to_noACR_riceACR_dic.keys()))
    for eachcate in store_broad_to_noACR_dic:
        cate_acr_dic = store_broad_to_noACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToNoACR' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_noACR_riceACR_num)
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##updating 110623
    ##conduct the enrichment test
    ##See if ACR under H3K27me3 in rice would like to be enriched on the ACR under the H3K27me3
    ##for the broad to broad
    store_final_line_list = []

    BtoB_H3K27rice_H3K27maize = len(list(store_broad_to_broad_riceACR_dic['H3K27me3'].keys()))
    BtoB_H3K27rice_absentmaize = len(list(store_broad_to_broad_riceACR_dic['H3K27me3Absent'].keys()))

    BtoNB_H3K27rice_H3K27maize = len(list(store_broad_to_notbroad_riceACR_dic['H3K27me3'].keys()))
    BtoNB_H3K27rice_absentmaize = len(list(store_broad_to_notbroad_riceACR_dic['H3K27me3Absent'].keys()))

    BtoNoACR_H3K27rice_H3K27maize = len(list(store_broad_to_noACR_dic['H3K27me3'].keys()))
    BtoNoACR_H3K27rice_absentmaize = len(list(store_broad_to_noACR_dic['H3K27me3Absent'].keys()))

    ##for the ACR in rice not overlapped with the H3K27me3
    BtoB_notH3K27rice_H3K27maize = len(list(store_broad_to_broad_riceACRnotCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoB_notH3K27rice_absentmaize = len(list(store_broad_to_broad_riceACRnotCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    BtoNB_notH3K27rice_H3K27maize = len(list(store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoNB_notH3K27rice_absentmaize = len(list(store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    BtoNoACR_notH3K27rice_H3K27maize = len(list(store_broad_to_noACR_notCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoNoACR_notH3K27rice_absentmaize = len(list(store_broad_to_noACR_notCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    ##for the notbroad rice ACR under H3K27me3 to all the other cases
    NBtoB_H3K27rice_H3K27maize = len(list(store_notbroad_to_broad_riceACR_dic['H3K27me3'].keys()))
    NBtoB_H3K27rice_absentmaize = len(list(store_notbroad_to_broad_riceACR_dic['H3K27me3Absent'].keys()))

    NBtoNB_H3K27rice_H3K27maize = len(list(store_notbroad_to_notbroad_riceACR_dic['H3K27me3'].keys()))
    NBtoNB_H3K27rice_absentmaize = len(list(store_notbroad_to_notbroad_riceACR_dic['H3K27me3Absent'].keys()))

    NBtoNoACR_H3K27rice_H3K27maize = len(list(store_notbroad_to_noACR_dic['H3K27me3'].keys()))
    NBtoNoACR_H3K27rice_absentmaize = len(list(store_notbroad_to_noACR_dic['H3K27me3Absent'].keys()))


    ##use the fisher to check the enrichment
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoB_H3K27rice_absentmaize],
                                            [BtoB_notH3K27rice_H3K27maize, BtoB_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoB' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoB_H3K27rice_absentmaize) + '\t' + str(BtoB_notH3K27rice_H3K27maize) + '\t' + str(BtoB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNB_H3K27rice_H3K27maize, BtoNB_H3K27rice_absentmaize],
                                            [BtoNB_notH3K27rice_H3K27maize, BtoNB_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNB' + '\t' + str(pvalue) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(
        BtoNB_H3K27rice_absentmaize) + '\t' + str(BtoNB_notH3K27rice_H3K27maize) + '\t' + str(
        BtoNB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNoACR_H3K27rice_H3K27maize, BtoNoACR_H3K27rice_absentmaize],
                                            [BtoNoACR_notH3K27rice_H3K27maize, BtoNoACR_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNoACR' + '\t' + str(pvalue) + '\t' + str(BtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_H3K27rice_absentmaize) + '\t' + str(BtoNoACR_notH3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_CompareTonotH3K27rice.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    store_final_line_list = []
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoB_H3K27rice_absentmaize],
                                            [NBtoB_H3K27rice_H3K27maize, NBtoB_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoB' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoB_H3K27rice_absentmaize) + '\t' + str(NBtoB_H3K27rice_H3K27maize) + '\t' + str(NBtoB_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNB_H3K27rice_H3K27maize, BtoNB_H3K27rice_absentmaize],
                                            [NBtoNB_H3K27rice_H3K27maize, NBtoNB_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNB' + '\t' + str(pvalue) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(
        BtoNB_H3K27rice_absentmaize) + '\t' + str(NBtoNB_H3K27rice_H3K27maize) + '\t' + str(
        NBtoNB_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNoACR_H3K27rice_H3K27maize, BtoNoACR_H3K27rice_absentmaize],
                                            [NBtoNoACR_H3K27rice_H3K27maize, NBtoNoACR_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNoACR' + '\t' + str(pvalue) + '\t' + str(BtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_H3K27rice_absentmaize) + '\t' + str(NBtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        NBtoNoACR_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_CompareToNotBroadH3K27me3rice.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##check if H3K27me3 covered BtoB
    store_final_line_list = []
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoNB_H3K27rice_H3K27maize],
                                            [BtoB_notH3K27rice_absentmaize, BtoNB_notH3K27rice_absentmaize]],
                                           alternative='less')

    final_line = 'H3K27compareNotH3K27me3' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(BtoB_notH3K27rice_absentmaize) + '\t' + \
        str(BtoNB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_H3K27compareNotH3K27me3.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##use the fisher to check the enrichment of broad compare to other and restricted ones










