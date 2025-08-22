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


def subfunction_summarize_overview (ipt_final_summary_fl, input_output_dir,spe1_prefix,spe2_prefix):


    ##eight groups
    store_shared_ACR_dic = {}
    store_var_ACR_dic = {}
    store_spes_ACR_dic = {}

    count = 0
    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1

            if count != 1:

                spe1_ACR = col[0]
                spe2_region = col[2]
                spe1_celltype = col[1]
                spe2_ACR = col[3]
                spe2_celltype = col[4]
                syn_region_id = col[5]

                if syn_region_id != 'none':

                    if spe2_region != 'none':

                        ##for the shared ACR
                        if spe2_ACR != 'none':

                            ACR_cate = 'na'

                            if spe1_celltype == 'broadly_accessible' and spe2_celltype == 'broadly_accessible':
                                ACR_cate = 'bACRbACR'
                            if spe1_celltype != 'broadly_accessible' and spe2_celltype != 'broadly_accessible':
                                ACR_cate = 'ctACRctACR'
                            if spe1_celltype == 'broadly_accessible' and spe2_celltype != 'broadly_accessible':
                                ACR_cate = 'bACRctACR'
                            if spe1_celltype != 'broadly_accessible' and spe2_celltype == 'broadly_accessible':
                                ACR_cate = 'ctACRbACR'

                            store_shared_ACR_dic[spe1_ACR] = ACR_cate

                        ##for the variable ACR
                        else:

                            if spe1_celltype == 'broadly_accessible':
                                ACR_cate = 'bACR'
                            else:
                                ACR_cate = 'ctACR'

                            store_var_ACR_dic[spe1_ACR] = ACR_cate

                    else:

                        ##for the species specific ACR
                        if spe1_celltype == 'broadly_accessible':
                            ACR_cate = 'bACR'
                        else:
                            ACR_cate = 'ctACR'

                        store_spes_ACR_dic[spe1_ACR] = ACR_cate

    store_temp_variable_ACR_dic = {}
    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:
                spe1_ACR = col[0]
                syn_region_id = col[5]

                if syn_region_id != 'none':

                    if ',' in syn_region_id:

                        final_cate = 'na' + '\t' + 'na'

                    else:
                        final_cate = 'wrong'
                        if spe1_ACR in store_shared_ACR_dic:
                            cate = store_shared_ACR_dic[spe1_ACR]
                            final_cate = 'shared_ACR' + '\t' +  cate
                        if spe1_ACR in store_var_ACR_dic:
                            cate = store_var_ACR_dic[spe1_ACR]
                            final_cate = 'variable_ACR' + '\t' + cate
                        if spe1_ACR in store_spes_ACR_dic:
                            cate = store_spes_ACR_dic[spe1_ACR]
                            final_cate = 'speciesspec_ACR' + '\t' + cate

                    final_line = eachline + '\t' + final_cate
                    store_final_line_list.append(final_line)

                else:
                    final_cate = 'na' + '\t' + 'na'
                    final_line = eachline + '\t' + final_cate
                    store_final_line_list.append(final_line)

                ##updating 082225
                all_col = final_line.strip().split()
                ACR_cate_infor = all_col[-2]
                if ACR_cate_infor == 'variable_ACR':
                    blasted_region_ACRinfor = col[2] + '__' + col[3]
                    if spe1_ACR in store_temp_variable_ACR_dic:
                        store_temp_variable_ACR_dic[spe1_ACR][blasted_region_ACRinfor] = 1
                    else:
                        store_temp_variable_ACR_dic[spe1_ACR] = {}
                        store_temp_variable_ACR_dic[spe1_ACR][blasted_region_ACRinfor] = 1


            else:
                first_line = eachline + '\t' + spe1_prefix + '_ACRcate1' + '\t' + spe1_prefix + '_ACRcate2'
                store_final_line_list.append(first_line)


    ##updating 082225
    ##delete the varialble ACR line if the ACR already has blasted to ACR
    store_target_spe1ACR_blastedregion_dic = {}
    for eachspe1_ACR in store_temp_variable_ACR_dic:
        blasted_region_ACRinfor_dic = store_temp_variable_ACR_dic[eachspe1_ACR]


        blasted_count = 0
        for eachregionACR in blasted_region_ACRinfor_dic:

            mt = re.match('(.+)__(.+)',eachregionACR)
            ACRinfor = mt.group(2)
            blastedregion = mt.group(1)

            if ACRinfor != 'none':
                blasted_count += 1

                target_spe1_ACR_blastedregion = eachspe1_ACR + '__' + blastedregion
                store_target_spe1ACR_blastedregion_dic[target_spe1_ACR_blastedregion] = 1

        if blasted_count == 0:

            for eachregionACR in blasted_region_ACRinfor_dic:
                mt = re.match('(.+)__(.+)', eachregionACR)
                blastedregion = mt.group(1)
                target_spe1_ACR_blastedregion = eachspe1_ACR + '__' + blastedregion

                store_target_spe1ACR_blastedregion_dic[target_spe1_ACR_blastedregion] = 1

    store_final_line_final_list = []
    for eachline in store_final_line_list:
        col = eachline.strip().split('\t')
        spe1ACR_blastedregion = col[0] + '__' + col[2]
        ACRcate = col[6]
        if ACRcate != 'variable_ACR':
            store_final_line_final_list.append(eachline)
        else:
            if spe1ACR_blastedregion in store_target_spe1ACR_blastedregion_dic:

                blastACR = col[3]
                if blastACR == 'none':
                    store_final_line_final_list.append(eachline)
                else:
                    final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\t' + 'shared_ACR' + '\t' + col[7]
                    store_final_line_final_list.append(final_line)

    with open (input_output_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '.syntenic.txt','w+') as opt:
        for eachline in store_final_line_final_list:
            opt.write(eachline + '\n')

    ##plot the comparison















