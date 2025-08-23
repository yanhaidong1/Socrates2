#!/usr/bin/env python

import re
import glob
import sys
import subprocess
import os

##updating 053124 we will make a Os.syntenic_genes.all_os_ACRs.bed file



def s1_subfunction_build_ACR_in_syn_region (spe1_ACR_fl,spe2_ACR_fl,input_syntenic_region_fl,opt_dir,
                                         target_spe1,target_spe2):


    store_target_synR_spe1_gene_dic = {}
    store_target_synR_spe2_gene_dic = {}
    count = 0
    with open (input_syntenic_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            count += 1
            if count != 1:

                syntenic_region_ID = col[-2]
                spe_nm = col[-3]

                spe1_syntenic_gene = col[3]
                spe2_syntenic_gene = col[7]

                spe1_chr = col[0]
                spe1_st = col[1]
                spe1_ed = col[2]

                spe2_chr = col[4]
                spe2_st = col[5]
                spe2_ed = col[6]

                if syntenic_region_ID in store_target_synR_spe1_gene_dic:
                    store_target_synR_spe1_gene_dic[syntenic_region_ID][spe1_syntenic_gene] = {'chr':spe1_chr,
                                                                                               'st':spe1_st,
                                                                                               'ed':spe1_ed}
                else:
                    store_target_synR_spe1_gene_dic[syntenic_region_ID] = {}
                    store_target_synR_spe1_gene_dic[syntenic_region_ID][spe1_syntenic_gene] = {'chr': spe1_chr,
                                                                                               'st': spe1_st,
                                                                                               'ed': spe1_ed}

                if spe_nm == target_spe2:

                    if syntenic_region_ID in store_target_synR_spe2_gene_dic:
                        store_target_synR_spe2_gene_dic[syntenic_region_ID][spe2_syntenic_gene] = {'chr':spe2_chr,
                                                                                                   'st':spe2_st,
                                                                                                   'ed':spe2_ed}
                    else:
                        store_target_synR_spe2_gene_dic[syntenic_region_ID] = {}
                        store_target_synR_spe2_gene_dic[syntenic_region_ID][spe2_syntenic_gene] = {'chr': spe2_chr,
                                                                                                   'st': spe2_st,
                                                                                                   'ed': spe2_ed}

    #print(store_target_synR_spe2_gene_dic)

    ##intersect with region for the spe1
    store_final_line_list = []
    for eachsyntenicID in store_target_synR_spe1_gene_dic:

        chr = ''
        all_loc_list = []
        for eachgene in store_target_synR_spe1_gene_dic[eachsyntenicID]:
            spe1_syntenic_gene_chr = store_target_synR_spe1_gene_dic[eachsyntenicID][eachgene]['chr']
            spe1_syntenic_gene_st = store_target_synR_spe1_gene_dic[eachsyntenicID][eachgene]['st']
            spe1_syntenic_gene_ed = store_target_synR_spe1_gene_dic[eachsyntenicID][eachgene]['ed']
            chr = spe1_syntenic_gene_chr

            all_loc_list.append(int(float(spe1_syntenic_gene_st)))
            all_loc_list.append(int(float(spe1_syntenic_gene_ed)))

        max_loc = max(all_loc_list)
        min_loc = min(all_loc_list)

        final_line = chr + '\t' + str(min_loc) + '\t' + str(max_loc) + '\t' + eachsyntenicID
        store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_' + target_spe1 + '_syntenic_region.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + target_spe1 + '_syntenic_region.txt > ' + \
          opt_dir + '/temp_' + target_spe1 + '_syntenic_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + target_spe1 + '_syntenic_region_sorted.txt' + \
          ' -b ' + spe1_ACR_fl + ' > ' + opt_dir + '/opt_intersect_' + target_spe1 + '_syn_ACR.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)


    print(store_target_synR_spe2_gene_dic)
    ##intersect with region for the spe2
    store_final_line_list = []
    for eachsyntenicID in store_target_synR_spe2_gene_dic:

        chr = ''
        all_loc_list = []
        for eachgene in store_target_synR_spe2_gene_dic[eachsyntenicID]:
            spe2_syntenic_gene_chr = store_target_synR_spe2_gene_dic[eachsyntenicID][eachgene]['chr']
            spe2_syntenic_gene_st = store_target_synR_spe2_gene_dic[eachsyntenicID][eachgene]['st']
            spe2_syntenic_gene_ed = store_target_synR_spe2_gene_dic[eachsyntenicID][eachgene]['ed']
            chr = spe2_syntenic_gene_chr

            all_loc_list.append(int(float(spe2_syntenic_gene_st)))
            all_loc_list.append(int(float(spe2_syntenic_gene_ed)))


        ##Here we consider all the gene body information

        max_loc = max(all_loc_list)
        min_loc = min(all_loc_list)

        final_line = chr + '\t' + str(min_loc) + '\t' + str(max_loc) + '\t' + eachsyntenicID
        store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_' + target_spe2 + '_syntenic_region.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + target_spe2 + '_syntenic_region.txt > ' + \
          opt_dir + '/temp_' + target_spe2 + '_syntenic_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + target_spe2 + '_syntenic_region_sorted.txt' + \
          ' -b ' + spe2_ACR_fl + ' > ' + opt_dir + '/opt_intersect_' + target_spe2 + '_syn_ACR.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)


def s1_subfunction_build_Os_syntenic_genes_all_os_ACR_fl (input_syntenic_region_fl,opt_dir,
                                                          target_spe1):


    store_acr_syntenic_ID_dic = {}
    with open (opt_dir + '/opt_intersect_' + target_spe1 + '_syn_ACR.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr_syntenic_region_ID = col[3]
            store_acr_syntenic_ID_dic[acr_syntenic_region_ID] = 1

    store_final_line_list = []
    with open (input_syntenic_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            synID = col[-2]

            if synID in store_acr_syntenic_ID_dic:
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + synID
                store_final_line_list.append(final_line)

    with open (opt_dir + '/' +target_spe1 + '.syntenic_genes.all_' + target_spe1 + '_ACRs.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

































