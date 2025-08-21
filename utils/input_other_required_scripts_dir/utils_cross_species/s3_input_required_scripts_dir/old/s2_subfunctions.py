#!/usr/bin/env python

##updating 112523 we will check the broad and CT peaks overlapping with 128 cell types
##updating 111723 we will add a new function to check the accessiblility of the rice altas peaks for the syntenic peaks
##updating 111523 we will add the categories of ACR nearby genes
##updating 111523 check if broad ACRs were relative higher enriched compared to the others


import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats


def subfunction_summarize_overview (ipt_final_summary_blast_fl,ipt_sep2_syntenic_region_fl,input_rice_phy_score_fl,opt_dir,
                                    ipt_spe1_ACR_gene_Cate_fl,ipt_spe2_ACR_gene_Cate_fl,ipt_spe1_ACR_celltype_fl,
                                    spe1_prefix,spe2_prefix,
                                    s2_s1_decide_celltype_cate,
                                    s2_s1_open_check_phy_score):


    store_syntenic_region_allRiceACR_dic = {}
    store_syntenic_region_RiceACR_MaizeACR_dic = {}
    store_allRiceACR_dic = {}
    store_allRiceACR_BlastToMaizeACR_dic = {}
    store_allRiceACR_BlastToMaizeRegionNotACR_dic = {}

    store_allMaizeACR_couldBeBlastedToRiceACR_dic = {}

    ##store CT acr and broad ACR
    store_allRiceACR_diffcate_dic = {}
    store_allRiceACR_BlastToMaizeACR_diffcate_dic = {}
    store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic = {}

    ##store maize acr
    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic = {}

    ##store rice acr
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic = {}

    ##store rice and maize acr at same time
    ##udpating 111523
    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic = {}


    ##store conserved broad and also the CT specific acr to check the score
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic = {}

    ##store broad ACR to maize different cell type
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic = {}


    ##store_riceACR_ct_cate_dic
    ##updating 112523
    store_riceACR_CT_cate_dic = {}

    count = 0
    with open (ipt_final_summary_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRID = col[0]
                riceACRloc = col[1]
                maizeACRloc = col[5]
                ricemaizeACRloc = riceACRloc + ',' + maizeACRloc

                acrline = '\t'.join(col[1].split('_'))
                store_riceACR_CT_cate_dic[acrline] = col[8]


                if s2_s1_decide_celltype_cate == 'CTVer':
                    rice_celltype = col[8]
                    maize_celltype = col[14]
                else:
                    if s2_s1_decide_celltype_cate == 'CoverVer':
                        rice_celltype = col[15]
                        maize_celltype = col[17]
                    else:
                        print('No rice and maize cell type identfied, please check the decide celltype cate is right or not')
                        break


                SyntenicRegionID_str = col[11]
                SyntenicRegionID_list = SyntenicRegionID_str.split(',')
                for eachSyntenicRegionID in SyntenicRegionID_list:
                    store_syntenic_region_allRiceACR_dic[eachSyntenicRegionID] = 1
                    store_allRiceACR_dic[riceACRloc] = 1

                    ##add cell type information
                    if s2_s1_decide_celltype_cate == 'CTVer':
                        if rice_celltype == 'broadly_accessible':
                            final_cate = 'broadACR'
                        else:
                            final_cate = 'restrictedACR'
                    else:

                        if s2_s1_decide_celltype_cate == 'CoverVer':
                            rice_celltype_cover_cate = col[16]
                            if rice_celltype_cover_cate == 'broad':
                                final_cate = 'broadACR'
                            else:
                                if rice_celltype_cover_cate == 'restricted':
                                    final_cate = 'restrictedACR'
                                else:
                                    final_cate = 'otherACR'

                        else:
                            print('please check the decide cell type cate is right or not')
                            break

                    if final_cate in store_allRiceACR_diffcate_dic:
                        store_allRiceACR_diffcate_dic[final_cate][riceACRloc] = 1
                    else:
                        store_allRiceACR_diffcate_dic[final_cate] = {}
                        store_allRiceACR_diffcate_dic[final_cate][riceACRloc] = 1


                    if riceACRID != 'none':

                        if maizeACRloc != 'none':

                            store_allRiceACR_BlastToMaizeACR_dic[riceACRloc] = 1
                            store_allMaizeACR_couldBeBlastedToRiceACR_dic[maizeACRloc] = 1

                            ##add cell type information
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype == 'broadly_accessible':
                                    final_rice_cate = 'broadACR'
                                else:
                                    final_rice_cate = 'restrictedACR'
                            else:
                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    rice_celltype_cover_cate = col[16]
                                    if rice_celltype_cover_cate == 'broad':
                                        final_rice_cate = 'broadACR'
                                    else:
                                        if rice_celltype_cover_cate == 'restricted':
                                            final_rice_cate = 'restrictedACR'
                                        else:
                                            final_rice_cate = 'otherACR'
                                else:
                                    print('please check the decide cell type cate is right or not')
                                    break


                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_diffcate_dic:
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate][riceACRloc] = 1

                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if maize_celltype == 'broadly_accessible':
                                    final_maize_cate = 'broadACR'
                                else:
                                    final_maize_cate = 'restrictedACR'
                            else:
                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    maize_celltype_cover_cate = col[18]
                                    if maize_celltype_cover_cate == 'broad':
                                        final_maize_cate = 'broadACR'
                                    else:
                                        if maize_celltype_cover_cate == 'restricted':
                                            final_maize_cate = 'restrictedACR'
                                        else:
                                            final_maize_cate = 'otherACR'
                                else:
                                    print('please check the decide cell type cate is right or not')
                                    break


                            ##calculate the number of maize ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][maizeACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][maizeACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][
                                    final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][
                                    final_maize_cate][maizeACRloc] = 1

                            ##calculate the number of rice ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1

                            ##updating 111523
                            ##store teh rice maize ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[
                                    final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate][ricemaizeACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate][ricemaizeACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][ricemaizeACRloc] = 1


                            ##updating 110223
                            ##we will check the ACR in different CT
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype != 'broadly_accessible':

                                    rice_celltype_list = rice_celltype.split(',')
                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachrice_celltype in rice_celltype_list:
                                        for eachmaize_celltype in maize_celltype_list:

                                            if eachrice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                                                if eachmaize_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype]:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype][riceACRloc] = 1
                                                else:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype] = {}
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                ##updating 110223
                                if rice_celltype == 'broadly_accessible':

                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachmaize_celltype in maize_celltype_list:

                                        if rice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic:
                                            if eachmaize_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[rice_celltype]:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype][riceACRloc] = 1
                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype][riceACRloc] = 1
                                        else:
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype] = {}
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype][eachmaize_celltype] = {}
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype][eachmaize_celltype][riceACRloc] = 1


                            if s2_s1_decide_celltype_cate == 'CoverVer':

                                rice_celltype_cover_cate = col[16]
                                if rice_celltype_cover_cate == 'restricted':

                                    rice_celltype_list = rice_celltype.split(',')
                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachrice_celltype in rice_celltype_list:
                                        for eachmaize_celltype in maize_celltype_list:

                                            if eachrice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                                                if eachmaize_celltype in \
                                                        store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                            eachrice_celltype]:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1
                                                else:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype] = {}
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype][eachmaize_celltype][riceACRloc] = 1



                            ##for the syntenic check
                            if eachSyntenicRegionID in store_syntenic_region_RiceACR_MaizeACR_dic:

                                if riceACRID in store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID]:
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1
                                else:
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID] = {}
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1

                            else:
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID] = {}
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID] = {}
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1



                        else:
                            store_allRiceACR_BlastToMaizeRegionNotACR_dic[riceACRloc] = 1

                            ##add cell type information
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype == 'broadly_accessible':
                                    final_cate = 'broadACR'
                                else:
                                    final_cate = 'restrictedACR'
                            else:

                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    rice_celltype_cover_cate = col[16]
                                    if rice_celltype_cover_cate == 'broad':
                                        final_cate = 'broadACR'
                                    else:
                                        if rice_celltype_cover_cate == 'restricted':
                                            final_cate = 'restrictedACR'
                                        else:
                                            final_cate = 'otherACR'


                            if final_cate in store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic:
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate] = {}
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate][riceACRloc] = 1


    store_total_syntenic_region_dic = {}
    with open (ipt_sep2_syntenic_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_total_syntenic_region_dic[col[3]] = 1

    syntenic_region_num = len(list(store_total_syntenic_region_dic.keys()))

    ##Here we will check if the syntenic region contains the acr
    syntenic_region_capturing_riceACR_dic = {}
    for eachregion in store_total_syntenic_region_dic:
        if eachregion in store_syntenic_region_allRiceACR_dic:
            syntenic_region_capturing_riceACR_dic[eachregion] = 1
    syntenic_region_real_num_capturing_riceACR = len(list(syntenic_region_capturing_riceACR_dic.keys()))

    ##check the diff between the real capture rice ACR
    for eachregion in store_syntenic_region_allRiceACR_dic:
        if eachregion not in syntenic_region_capturing_riceACR_dic:
            print(eachregion + ' has something wrong')

    ######
    ##opt1: check number of ACR and region for different case
    store_final_line_list = []
    syntenic_region_num_capturing_riceACR = len(list(store_syntenic_region_allRiceACR_dic.keys()))
    riceACR_num_within_syntenic_region = len(list(store_allRiceACR_dic.keys()))
    riceACR_num_within_syntenic_region_blastTomaizeACR = len(list(store_allRiceACR_BlastToMaizeACR_dic.keys()))
    riceACR_num_within_syntenic_region_blastTomaizeregionNotACR = len(list(store_allRiceACR_BlastToMaizeRegionNotACR_dic.keys()))

    final_line = 'syntenic_region_num' + '\t' + str(syntenic_region_num)
    store_final_line_list.append(final_line)
    final_line = 'syntenic_region_num_capturing_riceACR' + '\t' + str(syntenic_region_num_capturing_riceACR)
    store_final_line_list.append(final_line)
    final_line = 'syntenic_region_real_num_capturing_riceACR' + '\t' + str(syntenic_region_real_num_capturing_riceACR)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region' + '\t' + str(riceACR_num_within_syntenic_region)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region_blastTomaizeACR' + '\t' + str(riceACR_num_within_syntenic_region_blastTomaizeACR)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region_blastTomaizeregionNotACR' + '\t' + str(riceACR_num_within_syntenic_region_blastTomaizeregionNotACR)
    store_final_line_list.append(final_line)

    ##we will store the cell type specific ones
    for eachcate in store_allRiceACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region' + '\t' + str(len(list(store_allRiceACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)
    for eachcate in store_allRiceACR_BlastToMaizeACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region_blastTomaizeACR' + '\t' + str(len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)
    for eachcate in store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region_blastTomaizeregionNotACR' + '\t' + str(len(list(store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_ACRnum.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##updating 111723 we will make a version to store all the ACRs
    store_final_acr_three_cate_line = []
    store_final_acr_three_cate_addCTcate_line = []
    for eachacr in store_allRiceACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line = acr_line + '\t' + 'All'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_BlastToMaizeACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line =  acr_line + '\t' + 'SharedAcc'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line = acr_line + '\t' + 'SharedInAcc'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_dic:
        if eachacr not in store_allRiceACR_BlastToMaizeACR_dic:
            if eachacr not in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
                acr_line = '\t'.join(eachacr.split('_'))
                final_line = acr_line + '\t' + 'NotShared'
                store_final_acr_three_cate_line.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt','w+') as opt:
        for eachline in store_final_acr_three_cate_line:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt > ' + \
          opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##store_riceACR_CT_cate_dic
    store_final_line_list = []
    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr_line = col[0] + '\t' + col[1] + '\t' + col[2]
            CTcate = store_riceACR_CT_cate_dic[acr_line]
            final_line = eachline + '\t' + CTcate
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate.txt > ' + \
          opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##updating 111723
    ##we will generate an opt to store the cell type information
    ##this step will add the three category to the final acr file
    store_final_line_list = []
    with open (ipt_spe1_ACR_celltype_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]

            if acrnm in store_allRiceACR_dic:
                if acrnm not in store_allRiceACR_BlastToMaizeACR_dic:
                    if acrnm not in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
                        final_cate = 'NotShared'
                    else:
                        final_cate = 'SharedInAcc'

                else:
                    final_cate = 'SharedAcc'

            else:
                final_cate = 'NotBlast'


            mt = re.match('.+;(.+)',col[3])
            celltype = mt.group(1)

            celltype_list = celltype.split(',')
            for eachcelltype in celltype_list:

                final_line = col[0] + '_' + col[1] + '_' + col[2] + '\t' + eachcelltype + '\t' + final_cate
                store_final_line_list.append(final_line)

    with open ( opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_celltypes.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    ##updating 111523
    ##Here we will check if the broad syntenic broad were significantly enriched compared to the others
    ##use the fisher
    ##ACRtoACR_broad,  ACRtoACR_restricted
    ##ACRtoOthers_broad, ACRtoOthers_restricted

    ##oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])

    ACRtoACR_broad = len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic['broadACR'].keys()))
    ACRtoACR_restricted = len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic['restrictedACR'].keys()))

    ACRtoOthers_broad = len(list(store_allRiceACR_diffcate_dic['broadACR'].keys())) - ACRtoACR_broad
    ACRtoOthers_restricted = len(list(store_allRiceACR_diffcate_dic['restrictedACR'].keys())) - ACRtoACR_restricted

    oddsratio, pvalue = stats.fisher_exact([[ACRtoACR_broad, ACRtoACR_restricted], [ACRtoOthers_broad, ACRtoOthers_restricted]],alternative='greater')

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_ACRnum_check_ACRtoACR_enrich.txt','w+') as opt:
        final_line = 'pval' + '\t' + 'ACRtoACR_broad' + '\t' + 'ACRtoACR_restricted' + '\t' + 'ACRtoOthers_broad' + '\t' + 'ACRtoOthers_restricted'
        opt.write(final_line + '\n')
        final_line = str(pvalue) + '\t' + str(ACRtoACR_broad) + '\t' + str(ACRtoACR_restricted) + '\t' + str(ACRtoOthers_broad) + '\t' + str(ACRtoOthers_restricted)
        opt.write(final_line + '\n')










    ######
    ##opt2: check cate of ACR in maize to the cate of ACRs in rice
    ##updating 110123
    ##check how many restricted ACRs in rice would be the broad or restricted in maize
    ##updating 111523 we will add the category of close to genes
    store_spe1_acr_gene_cate_dic = {}
    with open (ipt_spe1_ACR_gene_Cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            cate = col[1]
            store_spe1_acr_gene_cate_dic[acrloc] = cate

    store_spe2_acr_gene_cate_dic = {}
    with open (ipt_spe2_ACR_gene_Cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            cate = col[1]
            store_spe2_acr_gene_cate_dic[acrloc] = cate


    store_final_line_list = []
    store_final_line_acr_ver_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe2_prefix + 'ACRnum'
    store_final_line_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[eachriceCTcate]:
            maize_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(maize_acr_dic.keys())))
            store_final_line_list.append(final_line)
            for eachmaizeacr in maize_acr_dic:
                final_line = '\t'.join(eachmaizeacr.split('_')) + '\t' + eachmaizeCTcate + '\t' + eachmaizeacr
                store_final_line_acr_ver_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRnum_celltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')
    with open (opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype.txt','w+') as opt:
        for eachline in store_final_line_acr_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype.txt > ' + \
          opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    store_final_line_acr_ver_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe1_prefix + 'ACRnum'
    store_final_line_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[eachriceCTcate]:
            rice_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(rice_acr_dic.keys())))
            store_final_line_list.append(final_line)
            for eachriceacr in rice_acr_dic:
                if eachriceacr in store_spe1_acr_gene_cate_dic:
                    riceacr_gene_cate = store_spe1_acr_gene_cate_dic[eachriceacr]
                else:
                    riceacr_gene_cate = 'none'
                final_line = '\t'.join(eachriceacr.split('_')) + '\t' + eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + riceacr_gene_cate
                store_final_line_acr_ver_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')
    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acr_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt > ' + \
          opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    ##udpating 111523
    ##opt2: check number and name of ACR pairs (for both Os and Zm)
    store_final_line_acrnum_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe1_prefix + 'ACRnum' + '\t' + spe2_prefix + 'ACRnum'
    store_final_line_acrnum_list.append(first_line)

    store_final_line_acrloc_ver_list = []
    first_line = spe1_prefix + 'ACRchr' + '\t' + spe1_prefix + 'ACRst' + '\t' + spe1_prefix + 'ACRed' + '\t' + \
                 spe1_prefix + 'CTcate' + '\t' + spe1_prefix + 'GeneCate' + '\t' + \
                 spe2_prefix + 'ACRchr' + '\t' + spe2_prefix + 'ACRst' + '\t' + spe2_prefix + 'ACRed' + '\t' + \
                 spe2_prefix + 'CTcate' + '\t' + spe2_prefix + 'GeneCate'
    store_final_line_acrloc_ver_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[eachriceCTcate]:

            ##store the number
            rice_maize_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(rice_maize_acr_dic.keys())))
            store_final_line_acrnum_list.append(final_line)

            ##store the loc
            for eachricemaizeacr in rice_maize_acr_dic:

                rice_acr = eachricemaizeacr.split(',')[0]
                maize_acr = eachricemaizeacr.split(',')[1]

                if rice_acr in store_spe1_acr_gene_cate_dic:
                    riceacr_gene_cate = store_spe1_acr_gene_cate_dic[rice_acr]
                else:
                    riceacr_gene_cate = 'none'

                if maize_acr in store_spe2_acr_gene_cate_dic:
                    maizeacr_gene_cate = store_spe2_acr_gene_cate_dic[maize_acr]
                else:
                    maizeacr_gene_cate = 'none'

                final_line = '\t'.join(rice_acr.split('_')) + '\t' + eachriceCTcate + '\t' + riceacr_gene_cate + '\t' + \
                             '\t'.join(maize_acr.split('_')) + '\t' + eachmaizeCTcate + '\t' + maizeacr_gene_cate
                store_final_line_acrloc_ver_list.append(final_line)

    with open(opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acrnum_list:
            opt.write(eachline + '\n')
    with open(opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acrloc_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt > ' + \
          opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)



    ##opt2: check cell type of ACRs in both rice and maize
    store_final_line_allacr_list = []
    store_final_line_acrcount_list = []
    for eachricecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
        for eachmaizecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype]:
            riceacr_dic =  store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype][eachmaizecelltype]
            riceacr_count = len(list(riceacr_dic.keys()))
            final_line = eachricecelltype + '\t' + eachmaizecelltype + '\t' + str(riceacr_count)
            store_final_line_acrcount_list.append(final_line)

            for eachacr in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype][eachmaizecelltype]:
                final_line = '\t'.join(eachacr.split('_')) + '\t' + eachricecelltype + '\t' + eachmaizecelltype
                store_final_line_allacr_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt','w+') as opt:
        for eachline in store_final_line_allacr_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_EachCelltype.txt','w+') as opt:
        for eachline in store_final_line_acrcount_list:
            opt.write(eachline + '\n')

    ##updating 110223
    ##opt2: check the broad ACRs in both rice and maize
    store_final_line_allacr_list = []
    store_final_line_acrcount_list = []
    for eachricecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic:
        for eachmaizecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype]:
            riceacr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype][eachmaizecelltype]
            riceacr_count = len(list(riceacr_dic.keys()))
            final_line = eachricecelltype + '\t' + eachmaizecelltype + '\t' + str(riceacr_count)
            store_final_line_acrcount_list.append(final_line)

            for eachacr in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype][eachmaizecelltype]:
                final_line = '\t'.join(eachacr.split('_')) + '\t' + eachricecelltype + '\t' + eachmaizecelltype
                store_final_line_allacr_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_BroadACRtoEachCelltype.txt','w+') as opt:
        for eachline in store_final_line_allacr_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_BroadACRtoEachCelltype.txt','w+') as opt:
        for eachline in store_final_line_acrcount_list:
            opt.write(eachline + '\n')


    #####################################################
    ##check the conesrvation score for these cate of ACRs
    if s2_s1_open_check_phy_score == 'yes':

        store_physcore_dir = opt_dir + '/store_physcore_dir'
        if not os.path.exists(store_physcore_dir):
            os.makedirs(store_physcore_dir)

        ipt_spe1_acr_cate_fl = opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'

        cmd = 'bedtools intersect -wa -wb -a ' + ipt_spe1_acr_cate_fl + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_celltypecate.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ipt_spe1_acr_celltype_fl = opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt'

        cmd = 'bedtools intersect -wa -wb -a ' + ipt_spe1_acr_celltype_fl + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_percelltype.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##updating 110823
        ##compare the phylo score between CT and Broad peaks for alls
        store_rice_ct_broad_line_dic = {}
        count = 0
        with open(ipt_final_summary_blast_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    riceACRID = col[0]
                    riceACRloc = col[1]
                    maizeACRloc = col[5]
                    ricecelltypecate = col[8]
                    maizecelltypecate = col[14]

                    if ricecelltypecate == 'broadly_accessible':
                        final_line = '\t'.join(riceACRloc.split('_')) + '\t' + ricecelltypecate
                        store_rice_ct_broad_line_dic[final_line] = 1
                    else:
                        final_line = '\t'.join(riceACRloc.split('_')) + '\t' + 'celltypeSpecific'
                        store_rice_ct_broad_line_dic[final_line] = 1

        with open (store_physcore_dir + '/temp_rice_broad_CT_acr.txt','w+') as opt:
            for eachline in store_rice_ct_broad_line_dic:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_physcore_dir + '/temp_rice_broad_CT_acr.txt > ' + \
              store_physcore_dir + '/temp_rice_broad_CT_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'bedtools intersect -wa -wb -a ' + store_physcore_dir + '/temp_rice_broad_CT_acr_sorted.txt' + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_broadCTACR.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)




    ######
    ##opt3: check one ACR in rice may correspond to how many ACRs in maize
    store_allRiceACR_blastToMaizeACR_num_dic = {}
    store_allRiceACR_blastToMaizeACR_line_list = []
    for eachsyntenic_region in store_syntenic_region_RiceACR_MaizeACR_dic:
        riceACRID_maizeACRloc_dic = store_syntenic_region_RiceACR_MaizeACR_dic[eachsyntenic_region]

        for eachriceACRID in riceACRID_maizeACRloc_dic:
            maizeACRloc_dic = riceACRID_maizeACRloc_dic[eachriceACRID]
            maizeACRnum = str(len(list(maizeACRloc_dic.keys())))
            maizeACRnum_name = 'BlastTo_' + maizeACRnum + '_' + spe2_prefix + 'ACRs'

            if maizeACRnum_name in store_allRiceACR_blastToMaizeACR_num_dic:
                store_allRiceACR_blastToMaizeACR_num_dic[maizeACRnum_name] += 1
            else:
                store_allRiceACR_blastToMaizeACR_num_dic[maizeACRnum_name] = 1

            maizeACRloc_str = ','.join(list(maizeACRloc_dic.keys()))

            num_maize_acr = len(list(maizeACRloc_dic.keys()))

            ##check if it is from the same or different col
            store_chrnm_dic = {}
            for eachmaizeacr in maizeACRloc_dic:
                mt = re.match('(.+)_.+_.+',eachmaizeacr)
                chrnm = mt.group(1)
                store_chrnm_dic[chrnm] = 1

            final_line = eachsyntenic_region + '\t' + eachriceACRID  + '\t' + str(num_maize_acr) + '\t' + str(len(store_chrnm_dic.keys())) + '\t' + maizeACRloc_str
            store_allRiceACR_blastToMaizeACR_line_list.append(final_line)


    with open (opt_dir + '/opt3_num_' + spe1_prefix + '_blastTo_diffnumOf_' + spe2_prefix + '.txt','w+') as opt:
        for eachcate in store_allRiceACR_blastToMaizeACR_num_dic:
            num_ACR = store_allRiceACR_blastToMaizeACR_num_dic[eachcate]
            final_line = eachcate + '\t' + str(num_ACR)
            opt.write(final_line + '\n')

    with open (opt_dir + '/opt3_' + spe1_prefix + '_blastTo_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_allRiceACR_blastToMaizeACR_line_list:
            opt.write(eachline + '\n')


    ##updating 110323
    ##Here we will generate different files plotting to the browser
    store_regions_ACR_to_upload_To_browser_dir = opt_dir + '/store_regions_ACR_to_upload_To_browser_dir'
    if not os.path.exists(store_regions_ACR_to_upload_To_browser_dir):
        os.makedirs(store_regions_ACR_to_upload_To_browser_dir)

    #riceACR_num_within_syntenic_region = len(list(store_allRiceACR_dic.keys()))
    #riceACR_num_within_syntenic_region_blastTomaizeACR = len(list(store_allRiceACR_BlastToMaizeACR_dic.keys()))
    #riceACR_num_within_syntenic_region_blastTomaizeregionNotACR = len(list(store_allRiceACR_BlastToMaizeRegionNotACR_dic.keys()))

    ##for all the rice ACRs in the syntenic regions
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion.txt','w+') as opt:
        for eachline in store_allRiceACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion.txt > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##for ACR in syntenic region blastToMaize ACR
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR.txt','w+') as opt:
        for eachline in store_allRiceACR_BlastToMaizeACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR.txt' + ' > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##we will check the spe2 ACR that could be blasted to rice ACR
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_'+ spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs.txt','w+') as opt:
        for eachline in store_allMaizeACR_couldBeBlastedToRiceACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_'+ spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs.txt > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


##updating 111723
def subfunction_check_acc_in_rice_atlas_peaks (ipt_rice_atlas_peak_coverage_fl,ipt_total_rice_atlas_acr_fl,ipt_rice_leaf_syntenic_fl,ipt_rice_leaf_syntenic_addCTcate_fl,
                                               opt_dir,
                                               spe1_prefix,spe2_prefix,
                                               s2_s2_organct_coverage_cutoff):


    ##s2_s2_organct_coverage_cutoff would be 2.5 for the rice atlas

    store_all_rice_atlas_peak_dic = {}
    with open (ipt_total_rice_atlas_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_all_rice_atlas_peak_dic[acrnm] = 1


    ##we will check the H3K27m3 for each of organ one by one
    store_peak_organct_val_dic = {}
    store_total_organct_dic = {}
    with open(ipt_rice_atlas_peak_coverage_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            peaknm = col[1]
            organct = col[0]
            coverage_val = col[2]
            posneg_cate = col[3]
            if peaknm in store_peak_organct_val_dic:
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}
            else:
                store_peak_organct_val_dic[peaknm] = {}
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}

            store_total_organct_dic[organct] = 1


    store_final_line_list = []
    store_eachpeak_celltype_str_num_dic = {}
    for eachpeak in store_peak_organct_val_dic:

        ##we want to allow the peaks are in the final peaks we called
        if eachpeak in store_all_rice_atlas_peak_dic:

            store_celltype_list = []
            store_cover_val_list = []

            for eachorganct in store_peak_organct_val_dic[eachpeak]:
                coverage_val = store_peak_organct_val_dic[eachpeak][eachorganct]['coverage_val']
                posneg_cate = store_peak_organct_val_dic[eachpeak][eachorganct]['posneg_cate']

                if float(coverage_val) >= float(s2_s2_organct_coverage_cutoff):

                    ##updating 050423
                    if posneg_cate == 'Positive':
                        store_celltype_list.append(eachorganct)
                        store_cover_val_list.append(coverage_val)

            store_celltype_str = ','.join(store_celltype_list)
            store_cover_val_str = ','.join(store_cover_val_list)

            if store_celltype_str != '':
                celltypenum = len(store_celltype_list)
                final_line = 'CutoffCover:' + s2_s2_organct_coverage_cutoff + '\t' + eachpeak + '\t' + store_celltype_str + '\t' + store_cover_val_str + '\t' + str(
                    celltypenum)
                store_final_line_list.append(final_line)
                store_eachpeak_celltype_str_num_dic[eachpeak] = {'celltypestr': store_celltype_str,
                                                                 'celltypenum': str(celltypenum)}

    with open(opt_dir + '/temp_summary_DA_peak_string_filterCoverage.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##check the overlapping with the single rice called peaks
    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_091823/output_dir_103123/step01_species_compare_add_cate_dir/store_rice_to_maize_dir
    ##Here we will use the three categories files
    cmd = 'bedtools intersect -wa -wb -a ' + ipt_rice_leaf_syntenic_fl + ' -b '  + ipt_total_rice_atlas_acr_fl + ' > ' + \
          opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    with open (opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            atlas_acr = col[4] + '_' + col[5] + '_' + col[6]

            celltypenum = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypenum']

            final_line = eachline + '\t' + celltypenum
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    cmd = 'bedtools intersect -wa -wb -a ' + ipt_rice_leaf_syntenic_addCTcate_fl + ' -b '  + ipt_total_rice_atlas_acr_fl + ' > ' + \
          opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas_addCTcate.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    with open (opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas_addCTcate.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            atlas_acr = col[5] + '_' + col[6] + '_' + col[7]

            celltypenum = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypenum']

            final_line = eachline + '\t' + celltypenum
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype_addCTcate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')













