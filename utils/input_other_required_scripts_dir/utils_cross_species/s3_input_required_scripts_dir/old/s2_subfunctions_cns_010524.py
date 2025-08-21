#!/usr/bin/env python

##updating 120423 we will do the enrichment test for the restricted ACRs
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
import random
from statistics import mean


##updating 120723
##build the CNS intersection for the lineage specific ACRs
def subfunction_check_CNS_overlap_density (ipt_spe1_cns_gff_fl,ipt_final_summary_fl,opt_dir):

    store_all_target_loc_fl_dir = opt_dir + '/store_all_target_loc_fl_dir'
    if not os.path.exists(store_all_target_loc_fl_dir):
        os.makedirs(store_all_target_loc_fl_dir)

    store_intersect_cns_fl_dir = opt_dir + '/store_intersect_cns_fl_dir'
    if not os.path.exists(store_intersect_cns_fl_dir):
        os.makedirs(store_intersect_cns_fl_dir)

    ##check the number of col in the gff file
    number_col = 0
    with open(ipt_spe1_cns_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open (ipt_spe1_cns_gff_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open (ipt_spe1_cns_gff_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_cns.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_cns.txt > ' + \
          opt_dir + '/opt_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##we will consider several cate
    store_shared_A_spe1_broadorCT_loc_dic = {}
    store_shared_I_spe1_broadorCT_loc_dic = {}
    store_not_shared_spe1_broadorCT_loc_dic = {}

    store_shared_I_spe1_CTcelltype_loc_dic = {}
    store_not_shared_spe1_CTcelltype_loc_dic = {}
    store_shared_A_spe1_CTcelltype_loc_dic = {}

    count = 0

    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:

                ##store the different cate of acr
                spe1_acrID = col[0]
                spe1_acrloc = col[1]
                spe2_acrloc = col[5]

                spe1_celltype = col[8]

                ##for the shared-A
                if spe1_acrloc != 'none' and spe2_acrloc != 'none':

                    if spe1_celltype == 'broadly_accessible':
                        broadorCT = 'Broad'
                    else:
                        broadorCT = 'CT'

                    if broadorCT in store_shared_A_spe1_broadorCT_loc_dic:
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_A_spe1_CTcelltype_loc_dic:
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1


                ##for the shared-I
                if spe1_acrID != 'none' and spe2_acrloc == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_I_spe1_CTcelltype_loc_dic:
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'

                    if broadorCT in store_shared_I_spe1_broadorCT_loc_dic:
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1




                ##for the not shared
                if spe1_acrID == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_not_shared_spe1_CTcelltype_loc_dic:
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'


                    if broadorCT in store_not_shared_spe1_broadorCT_loc_dic:
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1


    ##save the loc information
    for eachcate in store_shared_A_spe1_broadorCT_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr.txt','w+') as opt:
            for eachloc in store_shared_A_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    for eachcate in store_shared_I_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_shared_I_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

    for eachcate in store_not_shared_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_not_shared_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)


    for eachcelltype in store_shared_I_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_shared_I_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

    for eachcelltype in store_not_shared_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_not_shared_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    for eachcelltype in store_shared_A_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_shared_A_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


    all_target_fl_list = glob.glob(store_all_target_loc_fl_dir + '/*sorted.txt')

    store_final_line_list = []
    for eachfl in all_target_fl_list:

        mt = re.match('.+/(.+)',eachfl)
        flnm = mt.group(1)

        mt = re.match('opt_(.+)_acr_sorted\.txt',flnm)
        catenm = mt.group(1)

        store_all_acr_loc_dic = {}
        with open (eachfl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                store_all_acr_loc_dic[acrloc] = 1

        cmd = 'bedtools intersect -wa -wb -a ' + eachfl + \
              ' -b ' + opt_dir + '/opt_cns_sorted.txt > ' + \
              opt_dir + '/opt_intersect_' + catenm + '_acr_cns.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##calculate the num of cns per ACR
        store_acr_cns_num_dic = {}
        with open (opt_dir + '/opt_intersect_' + catenm + '_acr_cns.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                cnsloc = col[3] + '_' + col[4] + '_' + col[5]

                if acrloc in store_acr_cns_num_dic:
                    store_acr_cns_num_dic[acrloc][cnsloc] = 1
                else:
                    store_acr_cns_num_dic[acrloc] = {}
                    store_acr_cns_num_dic[acrloc][cnsloc] = 1


        for eachacr in store_all_acr_loc_dic:

            if eachacr in store_acr_cns_num_dic:
                cnsloc_num = len(list(store_acr_cns_num_dic[eachacr].keys()))
            else:
                cnsloc_num = 0

            final_line = catenm + '\t' + eachacr + '\t' + str(cnsloc_num)
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_cate_acr_cns_num.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')





def subfunction_check_enrichment_for_celltype (ipt_final_cate_acr_cns_num_fl,opt_dir):


    store_sharedI_celltype_finalcatenum_dic = {}
    store_not_shared_celltype_finalcatenum_dic = {}
    store_sharedA_celltype_finalcatenum_dic = {}
    with open (ipt_final_cate_acr_cns_num_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            catenm = col[0]
            cnsloc_num = col[2]

            if cnsloc_num != '0':
                final_cate = 'CoverCNS'
            else:
                final_cate = 'NotCoverCNS'

            if '_Broad' not in catenm and '_CT' not in catenm:

                if '_unknown_cells' not in catenm:

                    ##focuse on the shared_I
                    if 'shared_I_' in catenm:

                        if catenm in store_sharedI_celltype_finalcatenum_dic:

                            if final_cate in store_sharedI_celltype_finalcatenum_dic[catenm]:
                                store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_sharedI_celltype_finalcatenum_dic[catenm] = {}
                            store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] = 1


                    if 'not_shared' in catenm:

                        if catenm in store_not_shared_celltype_finalcatenum_dic:

                            if final_cate in store_not_shared_celltype_finalcatenum_dic[catenm]:
                                store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_not_shared_celltype_finalcatenum_dic[catenm] = {}
                            store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] = 1

                    if 'shared_A_' in catenm:

                        if catenm in store_sharedA_celltype_finalcatenum_dic:

                            if final_cate in store_sharedA_celltype_finalcatenum_dic[catenm]:
                                store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_sharedA_celltype_finalcatenum_dic[catenm] = {}
                            store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] = 1






    ##for the sharedI
    store_fisher_test_final_line_list = []
    store_plot_prop_final_line_list = []
    for target_celltype in store_sharedI_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_sharedI_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_sharedI_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_sharedI_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:

                other_celltype_CoverCNS_num = store_sharedI_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_sharedI_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num


        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num, other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')


        final_line = 'shared_I' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                       str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)


        ##build the proprotion
        final_line = 'shared_I' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(target_celltype_CoverCNS_num/ (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'shared_I' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(other_celltype_coverCNS_total_num/(other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)

    ##for the not shared
    for target_celltype in store_not_shared_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_not_shared_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_not_shared_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_not_shared_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:
                other_celltype_CoverCNS_num = store_not_shared_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_not_shared_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num

        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num,
                                                 other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')

        final_line = 'not_shared' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(
            target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                     str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)

        ##build the proprotion
        final_line = 'not_shared' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(
            target_celltype_CoverCNS_num / (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'not_shared' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(other_celltype_coverCNS_total_num / (
                    other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)


    ##for the shared-A
    for target_celltype in store_sharedA_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_sharedA_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_sharedA_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_sharedA_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:
                other_celltype_CoverCNS_num = store_sharedA_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_sharedA_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num

        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num,
                                                 other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')

        final_line = 'shared_A' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(
            target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                     str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)

        ##build the proprotion
        final_line = 'shared_A' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(
            target_celltype_CoverCNS_num / (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'shared_A' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(
            other_celltype_coverCNS_total_num / (
                        other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)



    with open (opt_dir + '/opt_fisher_test_allcelltype_enrich_on_CNS.txt','w+') as opt:
        for eachline in store_fisher_test_final_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_plot_R_prop_allcelltype_enrich_on_CNS.txt','w+') as opt:
        for eachline in store_plot_prop_final_line_list:
            opt.write(eachline + '\n')


    ##for the CNS check the shared-A
    ##what happens to these ACRs



##upadting 121923
def subfunction_check_varACRotherSpe_inRiceRegion_in_Atlas (ipt_final_summary_fl,ipt_rice_atlas_fl,ipt_spe2_cns_gff_fl,input_rice_atlas_acr_coverage_fl,
                                                            opt_dir,spe1_prefix,spe2_prefix,s2_s2_organct_coverage_cutoff):


    ##save the cell type num for the spe2
    ##we will check the H3K27m3 for each of organ one by one
    store_all_rice_atlas_peak_dic = {}
    with open (ipt_rice_atlas_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_all_rice_atlas_peak_dic[acrnm] = 1


    store_peak_organct_val_dic = {}
    store_total_organct_dic = {}
    with open(input_rice_atlas_acr_coverage_fl, 'r') as ipt:
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



    ##check the number of col in the gff file
    number_col = 0
    with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_' + spe2_prefix + '_cns.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_cns.txt > ' + \
          opt_dir + '/opt_' + spe2_prefix + '_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)




    store_blastedRegion_dic_in_rice_dic = {}
    store_spe1_acr_dic = {}
    count = 0
    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                ##store the different cate of acr
                spe1_acrID = col[0]
                spe1_acrloc = col[1]
                spe2_acrloc = col[5]
                spe2_region = col[4]

                if spe2_region != 'none':
                    spe2_region_loc = '\t'.join(spe2_region.split('_'))
                    store_blastedRegion_dic_in_rice_dic[spe2_region_loc] = 1

                spe1_acr = '\t'.join(spe1_acrloc.split('_'))
                store_spe1_acr_dic[spe1_acr] = 1

    with open (opt_dir + '/temp_' + spe1_prefix + '_acr.txt','w+') as opt:
        for eachline in store_spe1_acr_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt > ' + \
          opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    with open (opt_dir + '/temp_' + spe2_prefix+ '_region.txt','w+') as opt:
        for eachline in store_blastedRegion_dic_in_rice_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_region.txt > ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##intersect with the Atlas
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt' + \
          ' -b ' + ipt_rice_atlas_fl + ' > ' + opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_acrregion_atlas_acr_celltype_dic = {}
    with open (opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            region_loc = col[0] + '_' + col[1] + '_' + col[2]
            rice_atlas_acr = col[3] + '_' + col[4] + '_' + col[5]
            celltype = col[7]

            if region_loc in store_acrregion_atlas_acr_celltype_dic:
                store_acrregion_atlas_acr_celltype_dic[region_loc][rice_atlas_acr] = celltype
            else:
                store_acrregion_atlas_acr_celltype_dic[region_loc] = {}
                store_acrregion_atlas_acr_celltype_dic[region_loc][rice_atlas_acr] = celltype


    ##opt_dir + '/opt_cns_sorted.txt'
    #intersect the cns with the spe1 ACR and spe2 region
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_cns_sorted.txt > ' + opt_dir + '/temp_' + spe1_prefix + '_intersect_acr_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe1_acr_capture_cns_dic = {}
    with open (opt_dir + '/temp_' + spe1_prefix + '_intersect_acr_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            store_spe1_acr_capture_cns_dic[acrloc] = 1


    ##intersect the spe2 cns with the sp2 region
    cmd = 'bedtools intersect -wa -wb -a ' +  opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_' + spe2_prefix + '_cns_sorted.txt > ' + \
          opt_dir + '/temp_intersect_' + spe2_prefix + '_region_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe2_region_contain_cns_dic = {}
    with open (opt_dir + '/temp_intersect_' + spe2_prefix + '_region_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            regionloc = col[0] + '_' + col[1] + '_' + col[2]
            store_spe2_region_contain_cns_dic[regionloc] = 1


    store_final_line_list = []
    with open (ipt_final_summary_fl,'r') as ipt:
        count = 0
        with open(ipt_final_summary_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    ##store the different cate of acr
                    spe2_region = col[4]
                    spe1_acr = col[1]

                    if spe1_acr in store_spe1_acr_capture_cns_dic:
                        spe1_cover_cns = 'yes'
                    else:
                        spe1_cover_cns = 'no'


                    if spe2_region in store_spe2_region_contain_cns_dic:
                        spe2_cover_cns = 'yes'
                    else:
                        spe2_cover_cns = 'no'

                    if spe2_region in store_acrregion_atlas_acr_celltype_dic:

                        rice_atlas_acr_dic = store_acrregion_atlas_acr_celltype_dic[spe2_region]

                        rice_atlas_acr_str = ';'.join(list(rice_atlas_acr_dic.keys()))

                        rice_atlas_celltye_list = []
                        for eachacr in rice_atlas_acr_dic:
                            rice_atlas_celltye_list.append(rice_atlas_acr_dic[eachacr])

                        rice_atlas_celltype_str = ';'.join(rice_atlas_celltye_list)

                    else:
                        rice_atlas_acr_str = 'none'
                        rice_atlas_celltype_str = 'none'

                    ##  store_eachpeak_celltype_str_num_dic[eachpeak] = {'celltypestr': store_celltype_str,
                     #                                            'celltypenum': str(celltypenum)}

                    ##check the cell type number
                    if rice_atlas_acr_str != 'none':


                        celltypenum_list = []
                        rice_atlas_acr_list = rice_atlas_acr_str.split(';')
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltypenum = int(store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypenum'])
                                celltypenum_list.append(celltypenum)

                        final_celltypenum = str(int(mean(celltypenum_list)))

                        ##check the number of organs
                        store_organ_dic = {}
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltype_string = store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypestr']
                                celltype_list = celltype_string.split(',')
                                for eachcelltype in celltype_list:
                                    mt = re.match('(.+)\.+',eachcelltype)
                                    organnm = mt.group(1)
                                    store_organ_dic[organnm] = 1

                        organnum = str(len(list(store_organ_dic.keys())))

                        ##updating 122423
                        ##if we find the cell type str all contains the leaf we will see there is no type cover
                        if len(list(store_organ_dic.keys())) == 1:

                            if list(store_organ_dic.keys())[0] == 'leaf':
                                organnum = 'none'
                                final_celltypenum = 'none'
                                rice_atlas_acr_str = 'none'


                    else:
                        organnum = 'none'
                        final_celltypenum = 'none'


                    final_line = eachline + '\t' + rice_atlas_acr_str + '\t' + rice_atlas_celltype_str + '\t' + spe1_cover_cns + '\t' + spe2_cover_cns + '\t' + organnum + '\t' + final_celltypenum
                    store_final_line_list.append(final_line)

                else:
                    final_line = eachline + '\t' + 'Atlas_ACR' + '\t' + 'Atlas_ACR_celltype' + '\t' + spe1_prefix + '_ACR_capture_CNS' + '\t' + spe2_prefix + '_region_capture_CNS' + '\t' + 'capture_atlas_organ_num' + '\t' + 'capture_altas_celltypenum'
                    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_add_' + spe2_prefix + '_atlas_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')












