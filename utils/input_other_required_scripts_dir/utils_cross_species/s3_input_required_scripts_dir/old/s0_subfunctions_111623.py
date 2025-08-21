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

def subfunction_check_syntenic_blocks (ipt_spe1_all_os_ACRs_fl, ipt_spe2_os_ACRs_fl,input_spe1_H3K27me3_fl,opt_dir,spe1_prefix,spe2_prefix):

    ##We aim to use the ipt_spe2_os_ACRs_fl to find the region
    store_spe2_region_dic = {}
    with open (ipt_spe2_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            spe2_syn_region = col[3]
            store_spe2_region_dic[spe2_syn_region] = 1

    store_final_line_list = []
    store_region_ID_gene_pair_sted_dic = {}
    with open (ipt_spe1_all_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region_ID = col[-1]
            if syntenic_region_ID in store_spe2_region_dic:
                store_final_line_list.append(eachline)

                gene_line = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3]

                if syntenic_region_ID in store_region_ID_gene_pair_sted_dic:
                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                    #store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)
                else:

                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID] = {}
                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                    #store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_genes.all_os_ACRs.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##check the distance
    store_bug_line_list = []
    store_final_line_list = []
    for eachsyntenic_region_ID in store_region_ID_gene_pair_sted_dic:

        gene_line_dic = store_region_ID_gene_pair_sted_dic[eachsyntenic_region_ID]
        if len(list(gene_line_dic.keys())) == 2:

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))
                gene_chr = gene_col[0]

            store_location_list.sort()

            loc_2nd = store_location_list[1]
            loc_3rd = store_location_list[2]
            loc_1st = store_location_list[0]
            loc_4th = store_location_list[3]

            ##updating 102623 we will use the 1st and 4th
            #distance = abs(loc_3rd - loc_2nd)
            distance = abs(loc_4th - loc_1st)
            final_line = gene_chr + '\t' + str(loc_1st) + '\t' + str(loc_4th) + '\t' + eachsyntenic_region_ID + '\t' + str(distance)
            store_final_line_list.append(final_line)

        else:
            store_bug_line_list.append(eachsyntenic_region_ID + '\t' + '\t'.join(list(gene_line_dic.keys())))

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_bug_line.txt', 'w+') as opt:
        for eachline in store_bug_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance.txt > ' + \
          opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##check the number of ACRs within the gene pair
    cmd = 'cut -f 1-3 ' + input_spe1_H3K27me3_fl + ' > ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt > ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance_sorted.txt' + \
          ' > ' + opt_dir + '/temp_intersect_acr_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_syntenic_region_spe1_ACR_dic = {}
    store_intersected_acr_dic = {}
    with open (opt_dir + '/temp_intersect_acr_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            ACRnm = col[0] + '_' + col[1] + '_' + col[2]
            syntenic_regionID = col[6]
            store_intersected_acr_dic[ACRnm] = 1

            if syntenic_regionID in store_syntenic_region_spe1_ACR_dic:
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID][ACRnm] = 1
            else:
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID] = {}
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID][ACRnm] = 1

    with open (opt_dir + '/opt_ACRs_in_syntenic_regions.txt','w+') as opt:
        for eachacr in store_intersected_acr_dic:
            opt.write(eachacr + '\n')

    store_final_line_list = []
    for eachregion in store_syntenic_region_spe1_ACR_dic:
        ACR_str = ','.join(list(store_syntenic_region_spe1_ACR_dic[eachregion].keys()))
        ACR_num = len(list(store_syntenic_region_spe1_ACR_dic[eachregion].keys()))
        final_line = eachregion + '\t' + ACR_str + '\t' + str(ACR_num)
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_contain_ACRs.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

##updating 102023
def subfunction_build_syntenic_bed (ipt_spe1_all_os_ACRs_fl,input_spe1_H3K27me3_fl,opt_dir,spe1_prefix):

    store_final_line_list = []
    store_region_ID_gene_pair_sted_dic = {}
    with open(ipt_spe1_all_os_ACRs_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region_ID = col[-1]

            store_final_line_list.append(eachline)

            gene_line = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3]

            if syntenic_region_ID in store_region_ID_gene_pair_sted_dic:
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                # store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)
            else:
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID] = {}
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                # store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)


    ##check the distance
    store_final_line_list = []
    for eachsyntenic_region_ID in store_region_ID_gene_pair_sted_dic:

        gene_line_dic = store_region_ID_gene_pair_sted_dic[eachsyntenic_region_ID]
        if len(list(gene_line_dic.keys())) == 2:

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))
                gene_chr = gene_col[0]

            store_location_list.sort()

            loc_2nd = store_location_list[1]
            loc_3rd = store_location_list[2]
            loc_1st = store_location_list[0]
            loc_4th = store_location_list[3]

            distance = abs(loc_3rd - loc_2nd)
            final_line = gene_chr + '\t' + str(loc_1st) + '\t' + str(
                loc_4th) + '\t' + eachsyntenic_region_ID + '\t' + str(distance)
            store_final_line_list.append(final_line)

        ##We do not consider this case, as there are three gene cases showing the genes overlapping with each other
        else:
            ##if we meet three or more
            ##we still need to find the location

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                gene_chr = gene_col[0]

                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))

            store_location_list.sort()

            total_site_len = len(store_location_list)
            gene_num = total_site_len/2

            print(store_location_list)
            print(gene_num)

            for i in range(int(gene_num) - 1):

                inter_st_posi_num = i*2 + 2
                inter_ed_posi_num = i*2 + 3
                print(inter_st_posi_num)
                print(inter_ed_posi_num)

                inter_loc_st = store_location_list[inter_st_posi_num - 1]
                inter_loc_ed = store_location_list[inter_ed_posi_num - 1]

                distance = abs(inter_loc_ed - inter_loc_st)

                ##
                #final_line = gene_chr + '\t' + str(inter_loc_st) + '\t' + str(inter_loc_ed) + '\t' + eachsyntenic_region_ID + '_' + str(i) + '\t' + str(distance)
                #store_final_line_list.append(final_line)

            ##1 2 3 4 5 6
            ##23 and 45

            ##12 34 56 78
            ##2 3 and 4 5 and 6 7


    with open(opt_dir + '/opt_all' + spe1_prefix +  '_syntenic_region_add_distance.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_all' + spe1_prefix +  '_syntenic_region_add_distance.txt > ' + \
          opt_dir + '/opt_all' + spe1_prefix + '_syntenic_region_add_distance_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##intersect with the syntenic regions
    cmd = 'cut -f 1-3 ' + input_spe1_H3K27me3_fl + ' > ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt' + ' -b ' + \
          opt_dir + '/opt_all' + spe1_prefix + '_syntenic_region_add_distance_sorted.txt > ' + \
          opt_dir + '/temp_intersect_all_' + spe1_prefix + '_ACR_syntenic.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_peak_dic = {}
    with open (opt_dir + '/temp_intersect_all_' + spe1_prefix + '_ACR_syntenic.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            peaknm = col[0] + '_' + col[1] + '_' + col[2]
            store_peak_dic[peaknm] = 1

    with open (opt_dir + '/opt_' + spe1_prefix + '_peak_in_syntenic_region.txt','w+') as opt:
        for eachpeak in store_peak_dic:
            opt.write(eachpeak + '\n')

