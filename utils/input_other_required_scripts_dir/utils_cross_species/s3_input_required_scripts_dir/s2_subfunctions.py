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

                spe1_ACR = col[1]
                spe2_region = col[2]
                spe1_celltype = col[4]
                spe2_ACR = col[3]
                spe2_celltype = col[-1]
                syn_region_id = col[7]

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

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:
                spe1_ACR = col[1]
                syn_region_id = col[7]

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

            else:
                first_line = eachline + '\t' + spe1_prefix + '_ACRcate1' + '\t' + spe1_prefix + '_ACRcate2'
                store_final_line_list.append(first_line)

    with open (input_output_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '.syntenic.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')















