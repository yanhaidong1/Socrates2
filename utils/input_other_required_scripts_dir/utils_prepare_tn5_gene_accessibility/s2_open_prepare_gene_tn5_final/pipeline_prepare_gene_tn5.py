#!/usr/bin/env python

##this script we will prepare the gene acc files
#import re
#import glob
#import sys
import subprocess
#import os
#from Bio import SeqIO

##1) prepare gene bed
##this will generate two target files: one is for the generate gene body acc and another is for the annotation

##2) prepare the genome size file
##3) prepare gene sparse

#input_all_required_pipeline_scripts_dir = sys.argv[1]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/07_1_prepare_gene_acc_files_060421/input_all_required_pipeline_scripts_dir

################
##for the step01
#input_gene_gff_fl = sys.argv[2]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/resources/MSU_rice_genome_version/MSU_r7.gene.gff
#input_genome_fai_fl = sys.argv[3]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/resources/MSU_rice_genome_version/MSU_r7.fa.fai
#target_nm = sys.argv[4] ##gene or mRNA since some gff contains the mRNA
#prefix = sys.argv[5] ##riceMSUr7 or others

################
##for the step02
#input_genome_fl = sys.argv[6]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/resources/MSU_rice_genome_version/MSU_r7.fa

################
##for the step03
#input_unique_tn5_bed_fl = sys.argv[7]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/02_keep_ref_mPing_tn5_comb_uniq_tn5_rmblack_051321/02_remove_black_list_052121/output_dir/opt_all_tn5_rmblack.bed

#input_output_dir = sys.argv[8]

def step01_prepare_gene_bed (input_other_required_scripts_dir,input_output_dir,
                             input_gene_gff_fl,input_genome_fai_fl,category):

    ##run the s1
    #store_s1_dir = input_output_dir + '/01_generate_gene_bed_fl_dir'
    #if not os.path.exists(store_s1_dir):
    #    os.makedirs(store_s1_dir)

    store_s1_dir = input_output_dir

    input_target_required_pipeline_scripts_dir = input_other_required_scripts_dir

    prepare_gene_bed_script = input_target_required_pipeline_scripts_dir + '/step01_prepare_gene_bed_fl.py'

    cmd = 'python ' + prepare_gene_bed_script + \
          ' ' + input_gene_gff_fl + \
          ' ' + input_genome_fai_fl + \
          ' ' + store_s1_dir + \
          ' ' + category
    print(cmd)
    subprocess.call(cmd,shell=True)

def step02_prepare_genome_size (input_other_required_scripts_dir,input_genome_fl,input_output_dir):

    ##run the s1
    #store_s1_dir = input_output_dir + '/02_generate_genome_size_fl_dir'
    #if not os.path.exists(store_s1_dir):
    #    os.makedirs(store_s1_dir)

    store_s1_dir = input_output_dir


    input_target_required_pipeline_scripts_dir = input_other_required_scripts_dir

    genome_size_script = input_target_required_pipeline_scripts_dir + '/step02_generate_genome_size_fl.py'

    cmd = 'python ' + genome_size_script + \
          ' ' + input_genome_fl + \
          ' ' + store_s1_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

def step03_prepare_gene_sparse (input_other_required_scripts_dir,input_tn5_bed_fl,input_genome_fai_fl,
                                input_output_dir):

    ##only step03 contain a gene sparse
    #step03_prepare_gene_sparse_dir = input_output_dir + '/step03_prepare_gene_sparse_dir'
    #if not os.path.exists(step03_prepare_gene_sparse_dir):
    #    os.makedirs(step03_prepare_gene_sparse_dir)


    ##s1 generate a nbinary sparse
    #store_s1_dir = step03_prepare_gene_sparse_dir + '/01_prepare_nonbinary_gene_sparse_fl_dir'
    #if not os.path.exists(store_s1_dir):
    #    os.makedirs(store_s1_dir)

    store_s1_dir = input_output_dir


    input_target_required_pipeline_scripts_dir = input_other_required_scripts_dir

    prepare_sparse_script = input_target_required_pipeline_scripts_dir + '/step03_prepare_nonbinary_gene_sparse_fl.py'
    nbfastSparsetn5_pl_script = input_target_required_pipeline_scripts_dir + '/fastSparse.nonbinary.peak.pl'

    input_gene_bed_fl = input_output_dir + '/opt_genes_500bpTSS_sorted.bed'

    cmd = 'python ' + prepare_sparse_script + \
          ' ' + input_tn5_bed_fl + \
          ' ' + input_gene_bed_fl + \
          ' ' + input_genome_fai_fl + \
          ' ' + store_s1_dir + \
          ' ' + nbfastSparsetn5_pl_script
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##s2 do the TPO normalization
    #store_s2_dir = step03_prepare_gene_sparse_dir + '/02_TPO_norm_dir'
    #if not os.path.exists(store_s2_dir):
    #    os.makedirs(store_s2_dir)

    #tpo_norm_script = input_target_required_pipeline_scripts_dir + '/step03_TPO_norm.R'
    #gene_sparse_fl = store_s1_dir + '/opt_gene_sparse_bed_dir/opt_gene_chnm.sparse'
    #gene_length_fl = input_output_dir + '/01_generate_gene_bed_fl_dir/opt_' + prefix + '_gene_length.txt'

    #cmd = 'Rscript ' + tpo_norm_script + \
    #      ' ' + gene_sparse_fl + \
    #      ' ' + gene_length_fl + \
    #      ' ' + store_s2_dir
    #subprocess.call(cmd,shell=True)


##updating 010425 save it to the object
def step04_save_to_soc_obj (input_other_required_scripts_dir,input_soc_obj_fl,input_gene_tn5_sparse_fl,input_gene_500bp_sorted_fl,input_prefix,input_output_dir):


    input_store_obj_R_script = input_other_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s2_open_prepare_gene_tn5_final/store_files_to_object.R'

    cmd = 'Rscript ' + input_store_obj_R_script + \
          ' ' + input_soc_obj_fl + \
          ' ' + input_gene_tn5_sparse_fl + \
          ' ' + input_gene_500bp_sorted_fl + \
          ' ' + input_prefix + \
          ' ' + input_output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)




def prepare_gene_tn5 (input_other_required_scripts_dir,input_tn5_bed_fl,input_output_dir,
                      input_gene_gff_fl,input_genome_fai_fl,
                      input_soc_obj_fl,input_prefix,
                      category = 'gene'):


    input_final_scripts_dir = input_other_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s2_open_prepare_gene_tn5_final'

    step01_prepare_gene_bed (input_final_scripts_dir,input_output_dir,
                                 input_gene_gff_fl,input_genome_fai_fl,category)

    #step02_prepare_genome_size (input_final_scripts_dir,input_genome_fl,input_output_dir)

    step03_prepare_gene_sparse (input_final_scripts_dir,input_tn5_bed_fl,input_genome_fai_fl,
                                input_output_dir)



    ##updating 010425
    input_gene_tn5_sparse_fl = input_output_dir + '/opt_gene_tn5.sparse'
    input_gene_500bp_sorted_fl = input_output_dir + '/opt_genes_500bpTSS_sorted.bed'
    step04_save_to_soc_obj (input_other_required_scripts_dir,input_soc_obj_fl,input_gene_tn5_sparse_fl,input_gene_500bp_sorted_fl,input_prefix,input_output_dir)


#prepare_gene_tn5 (input_other_required_scripts_dir,input_gene_gff_fl,input_genome_fai_fl,
#              target_nm,prefix,input_genome_fl,input_unique_tn5_bed_fl,input_output_dir)




