#!/usr/bin/env python

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
from Bio import SeqIO
import numpy as np



def multi_run(args):
    return blast_ACR_to_region (*args)

def blast_ACR_to_region(ipt_target_region_list,ipt_spe1_syn_region_ACR_fl,
                        ipt_spe2_syn_region_fl, ipt_spe1_genome_fa_fl,ipt_spe2_genome_fa_fl,
                        opt_dir):

    store_final_line_list = []
    for eachsyn_region_ID in ipt_target_region_list:

        store_target_ACR_loc_list = []
        with open (ipt_spe1_syn_region_ACR_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                ##the output will be chr + '\t' + str(min_loc) + '\t' + str(max_loc) + '\t' + \
                # eachsyntenicID + '\t' + 'ACRchr + '\t' + 'ACRst' + '\t' + 'ACRed'

                syntenicID = col[3]
                ACRchr = col[4]
                ACRst = col[5]
                ACRed = col[6]

                if syntenicID == eachsyn_region_ID:

                    final_line = ACRchr + '\t' + ACRst + '\t' + ACRed + '\t' + ACRchr + '_' + ACRst + '_' + ACRed
                    store_target_ACR_loc_list.append(final_line)

        with open (opt_dir + '/temp_ACR_loc.txt','w+') as opt:
            for eachline in store_target_ACR_loc_list:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_ACR_loc.txt > ' + \
              opt_dir + '/temp_ACR_loc_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##get the ACR fasta
        cmd = 'bedtools getfasta -fi ' + ipt_spe1_genome_fa_fl + \
              ' -bed ' + opt_dir + '/temp_ACR_loc_sorted.txt' + ' -name > ' + \
              opt_dir + '/temp_ACR.fa'
        print(cmd)
        subprocess.call(cmd,shell=True)


        ##chr + '\t' + str(min_loc) + '\t' + str(max_loc) + '\t' + eachsyntenicID
        store_target_spe2_syn_loc_list = []
        with open (ipt_spe2_syn_region_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                syntenicID = col[3]
                if syntenicID == eachsyn_region_ID:
                    store_target_spe2_syn_loc_list.append(eachline)

        with open (opt_dir + '/temp_syn_spe2_loc.txt','w+') as opt:
            for eachline in store_target_spe2_syn_loc_list:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_syn_spe2_loc.txt > ' + \
              opt_dir + '/temp_syn_spe2_loc_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##get the region fasta
        cmd = 'bedtools getfasta -fi ' + ipt_spe2_genome_fa_fl + \
              ' -bed ' + opt_dir + '/temp_syn_spe2_loc_sorted.txt' + ' -name > ' + \
              opt_dir + '/temp_syn_spe2.fa'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##make the blastdb
        ##udpating 072925 there is something wrong with parse seqids
        ##we will not use it
        #cmd = 'makeblastdb -in ' + opt_dir + '/temp_syn_spe2.fa' + ' -dbtype nucl -parse_seqids -out ' + opt_dir + '/temp_syn_spe2.db'
        #print(cmd)
        #subprocess.call(cmd,shell=True)
        cmd = 'makeblastdb -in ' + opt_dir + '/temp_syn_spe2.fa' + ' -dbtype nucl -out ' + opt_dir + '/temp_syn_spe2.db'
        print(cmd)
        subprocess.call(cmd,shell=True)


        ##blast
        ##blastn -query {input.species_acrs} -db {input.syntenic_region_db} -task blastn-short -out \
        #{output.blast} -evalue 1e-3 -max_target_seqs 4 -num_threads {threads} -word_size 7 -gapopen 5 -gapextend 2 \
        #-penalty -1 -reward 1 -dust no -outfmt 6

        cmd = 'blastn -query ' + opt_dir + '/temp_ACR.fa' + ' -db ' +  opt_dir + '/temp_syn_spe2.db' + \
              ' -task blastn-short -out ' + opt_dir + '/opt_blast_res.txt' + \
              ' -evalue 1e-3' + \
              ' -max_target_seqs 4' + \
              ' -num_threads 1' + \
              ' -word_size 7' + \
              ' -gapopen 5' + \
              ' -gapextend 2' + \
              ' -penalty -1 -reward 1 -dust no -outfmt 6'
        print(cmd)
        subprocess.call(cmd,shell=True)

        with open (opt_dir + '/opt_blast_res.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                store_final_line_list.append(eachline)

    with open (opt_dir + '/opt_final_blast_res.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

def subfunction_run_parallel_blast (opt_dir,ipt_all_region_fl,ipt_spe1_syn_region_ACR_fl,ipt_spe2_syn_region_fl,
                                    ipt_spe1_genome_fa_fl,ipt_spe2_genome_fa_fl,
                                    input_core_num):


    temp_all_output_dir = opt_dir + '/temp_all_output_dir'
    if not os.path.exists(temp_all_output_dir):
        os.makedirs(temp_all_output_dir)

    store_all_region_ID_list = []
    with open (ipt_all_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            store_all_region_ID_list.append(eachline)


    for x in range(0, int(input_core_num)):
        dir_code = x + 1
        temp_output_dir = temp_all_output_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    store_core_dic = {}
    core_count = -1
    array_split_list = np.array_split(store_all_region_ID_list, int(input_core_num))
    for eacharray in array_split_list:
        core_count += 1

        ##the array list contains cluster information
        array_list = []
        for eachitem in eacharray:
            array_list.append(eachitem)

        store_core_dic[str(core_count)] = array_list

    temp_all_output_dir_list = glob.glob(temp_all_output_dir + '/*')

    ##blast_ACR_to_region(ipt_target_region_list,ipt_spe1_syn_region_ACR_fl,
    #                    ipt_spe2_syn_region_fl, ipt_spe1_genome_fa_fl,ipt_spe2_genome_fa_fl,
    #                    opt_dir)


    pool = Pool(int(input_core_num))
    run_list = []
    for x in range(0, int(input_core_num)):
        each_func_argument = (store_core_dic[str(x)],
                              ipt_spe1_syn_region_ACR_fl,
                              ipt_spe2_syn_region_fl,
                              ipt_spe1_genome_fa_fl,
                              ipt_spe2_genome_fa_fl,
                              temp_all_output_dir_list[x])
        run_list.append(each_func_argument)
    pool.map(multi_run, run_list)



















