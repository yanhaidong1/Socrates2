#!/usr/bin/env python

##updating 020325 we will consider iD if there is no Name in the gff file
##updating 070821 consider the case that gene name has other name information change the geneID to the gene_feat
##updating 060421 prepare a file for the rice case and also prepare a gene length file for the step03
##updation 120920 generate a file that overlap 2kb upstream to 500-bp downstream
##updation 120320 add the chr besides the name
##this script we will prepare the gene annotation for the add_05 or other steps
import re
import glob
import sys
import subprocess
import os

input_gene_gff_fl = sys.argv[1]
input_black_chr_str = sys.argv[2]
input_genome_fai_fl = sys.argv[3]
input_output_dir = sys.argv[4]

target_nm = sys.argv[5] ##gene or mRNA since some gff contains the mRNA
#prefix = sys.argv[5] ##riceMSUr7 or others

def prepare_gene_bed_fl (input_gene_gff_fl,input_black_chr_str,input_genome_fai_fl,input_output_dir,target_nm):

    input_black_chr_str_list = input_black_chr_str.split(',')

    store_genome_name_dic = {}
    with open (input_genome_fai_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if col[0] not in input_black_chr_str_list:
                store_genome_name_dic[col[0]] = 1

    print('The following chromosomes are considering in the downstream analysis')
    print(store_genome_name_dic)

    store_gene_length_line_list = []
    store_gene_annotation_bed_list = []
    store_500TSS_bed_line_list = []
    with open (input_gene_gff_fl,'r') as ipt:
        for eachline in ipt:
            if not eachline.startswith('#'):
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                if col[2] == target_nm:

                    if col[0] in store_genome_name_dic:

                        start = col[3]
                        end = col[4]
                        dir = col[6]

                        if dir == '+':
                            real_end = end

                            if (int(start) - 500) < 0:
                                real_start = '0'
                            else:
                                real_start = str(int(start) - 500)
                        else:
                            real_start = start
                            real_end = str(int(end) + 500)

                        annot_col = col[8].split(';')

                        if len(annot_col) > 1:
                            if re.match('Name=(.+)',annot_col[1]):
                                mt = re.match('Name=(.+)',annot_col[1])
                                gene_feat = mt.group(1)
                            else:
                                mt = re.match('ID=(.+)', annot_col[0])
                                gene_feat = mt.group(1)
                        else:
                            mt = re.match('ID=(.+)', annot_col[0])
                            gene_feat = mt.group(1)



                        ##updating 070821 do not consider the  Pt and Mt
                        #if col[0] != 'ChrPt' and col[0] != 'ChrMt':
                        #    mt = re.match('Name=(.+)',annot_col[1])
                        #    gene_feat = mt.group(1)

                        #final_line = col[0] + '\t' + real_start + '\t' + real_end + '\t' + \
                        #             col[6] + '\t' + 'gene' + '\t' + gene_feat + '\t' + gene_id + '\t' + gene_feat
                        final_line = col[0] + '\t' + real_start + '\t' + real_end + '\t' + \
                                     col[6] + '\t' + 'gene' + '\t' + gene_feat + '\t' + gene_feat + '\t' + gene_feat

                        store_gene_annotation_bed_list.append(final_line)

                        final_line = col[0] + '\t' + real_start + '\t' + real_end + '\t' + \
                                     gene_feat + '\t' + col[6]
                        store_500TSS_bed_line_list.append(final_line)

                        ##updating prepare the gene length
                        length = int(end) - int(start)
                        #final_line = gene_id + '\t' + str(length)
                        final_line = gene_feat + '\t' + str(length)
                        store_gene_length_line_list.append(final_line)

    #with open (input_output_dir + '/opt_gene_length.txt','w+') as opt:
    #    for eachline in store_gene_length_line_list:
    #        opt.write(eachline + '\n')

    #with open (input_output_dir + '/opt_geneAnnotation.bed','w+') as opt:
    #    for eachline in store_gene_annotation_bed_list:
    #        opt.write(eachline + '\n')

    with open (input_output_dir + '/opt_genes_500bpTSS.bed','w+') as opt:
        for eachline in store_500TSS_bed_line_list:
            opt.write(eachline + '\n')

    ##we need to sort the bed information
    #cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/opt_geneAnnotation.bed' + ' > ' +  input_output_dir + '/opt_geneAnnotation_sorted.bed'
    #subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/opt_genes_500bpTSS.bed > ' + input_output_dir + '/opt_genes_500bpTSS_sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'rm ' + input_output_dir + '/opt_genes_500bpTSS.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##updation 120920
    #store_2000up500down_bed_line_list = []
    #with open (input_gene_gff_fl,'r') as ipt:
    #    for eachline in ipt:
    #        if not eachline.startswith('#'):
    #            eachline = eachline.strip('\n')
    #            col = eachline.strip().split()
    #            if col[2] == 'gene':

    #                if col[0] in store_genome_name_dic:

    #                    start = col[3]
    #                    end = col[4]
    #                    dir = col[6]

    #                    if dir == '+':
    #                        real_end = str(int(end) + 500)
    #                        if (int(start) - 2000) < 0:
    #                            real_start = '0'
    #                        else:
    #                            real_start = str(int(start) - 2000)
    #                    else:
    #                        if int(start) - 500 < 0:
    #                            real_start = '0'
    #                        else:
    #                            real_start = str(int(start) - 500)
    #                        real_end = str(int(end) + 2000)
                        
                        ##updating 070821 do not consider the  Pt and Mt
    #                    if col[0] != 'ChrPt' and col[0] != 'ChrMt':

    #                        annot_col = col[8].split(';')
    #                        mt = re.match('ID=(.+)',annot_col[0])
    #                        gene_id = mt.group(1)
    #                        mt = re.match('Name=(.+)',annot_col[1])
    #                        gene_feat = mt.group(1)

    #                        final_line = 'chr' + col[0] + '\t' + real_start + '\t' + real_end + '\t' + \
    #                                     gene_feat + '\t' + col[6]
    #                        store_2000up500down_bed_line_list.append(final_line)

    #with open (input_output_dir + '/opt_genes_2000up500down.bed','w+') as opt:
    #    for eachline in store_2000up500down_bed_line_list:
    #        opt.write(eachline + '\n')

    #cmd = 'sort -k1,1V -k2,2n ' + input_output_dir + '/opt_genes_2000up500down.bed > ' + input_output_dir + '/opt_genes_2000up500down_sorted.bed'
    #subprocess.call(cmd,shell=True)

prepare_gene_bed_fl (input_gene_gff_fl,input_black_chr_str,input_genome_fai_fl,input_output_dir,target_nm)
