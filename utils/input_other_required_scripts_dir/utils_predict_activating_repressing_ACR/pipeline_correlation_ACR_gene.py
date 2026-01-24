#!/usr/bin/env python

##this script is to perform the correlation between ACRs and genes
import re
import glob
import sys
import subprocess
import os


input_required_script_dir = sys.argv[1]

input_acr_fl = sys.argv[2]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/05_cross_species_072125/input_all_spe_celltype_acr_dir_updatePg/Pg_acr_celltype.txt

input_gene_gff_fl = sys.argv[3]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/resource/raw_data/Tifleaf3.gff3

input_scrna_object_fl = sys.argv[4]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/06_integrate_scRNAseq_072425/02_scrna_marker_enrichment_ver_073025/resource/seob_resct.updated.rds

##step03
input_acr_cpm_fl = sys.argv[5]
##/public2/home/yanhaidong/working_dir/yanhaidong/Pearlmillet_scATACseq/06_integrate_scRNAseq_072425/03_correlation_to_find_potential_enhancer_silencer_080525/resource/opt_perM_peaks_accessibility_cell_identity.txt

##step04
#input_atac_meta_fl = sys.argv[6]

input_select_pair_id_list_fl = sys.argv[6]
##if nothing we will use the na to represent

input_core_num = sys.argv[7]

input_output_dir = sys.argv[8]

###########
##add an additional one to include the ipt_target_colnm








###################
#final acr_pairs_fl
#Gm01    11311   11812   Gm01_11311_11812        Gm01    27344   28430   Glyma.01G000100.Wm82.a4.v1      .       -       15783
#Gm01    11881   12382   Gm01_11881_12382        Gm01    27344   28430   Glyma.01G000100.Wm82.a4.v1      .       -       15213
#Gm01    12495   12996   Gm01_12495_12996        Gm01    27344   28430   Glyma.01G000100.Wm82.a4.v1      .       -       14599
#Gm01    24658   25159   Gm01_24658_25159        Gm01    27344   28430   Glyma.01G000100.Wm82.a4.v1      .       -       2436
#Gm01    25962   26463   Gm01_25962_26463        Gm01    27344   28430   Glyma.01G000100.Wm82.a4.v1      .       -       1132

#####################
##final acr2gene file
##USE na if type1 and type2 do not have

#peak_id gene_id.close   dis2gene        peak.type1      peak.type2      peak2gene_id
#1       Gm01_11311_11812        Glyma.01G000100 15783   dACR    dACR    Gm01_11311_11812:Glyma.01G000100
#2       Gm01_11881_12382        Glyma.01G000100 15213   dACR    dACR    Gm01_11881_12382:Glyma.01G000100
#3       Gm01_12495_12996        Glyma.01G000100 14599   dACR    dACR    Gm01_12495_12996:Glyma.01G000100
#4       Gm01_24658_25159        Glyma.01G000100 2436    dACR    dACR    Gm01_24658_25159:Glyma.01G000100
#5       Gm01_25962_26463        Glyma.01G000100 1132    pACR    pACR    Gm01_25962_26463:Glyma.01G000100
#6       Gm01_29629_30130        Glyma.01G000100 -1450   pACR    pACR    Gm01_29629_30130:Glyma.01G000100
#7       Gm01_30174_30675        Glyma.01G000100 -1995   pACR    pACR    Gm01_30174_30675:Glyma.01G000100
#8       Gm01_46414_46915        Glyma.01G000137 5821    dACR    dACR    Gm01_46414_46915:Glyma.01G000137
#9       Gm01_49302_49803        Glyma.01G000137 2933    dACR    dACR    Gm01_49302_49803:Glyma.01G000137
#10      Gm01_78350_78851        Glyma.01G000322 -874    pACR    pACR    Gm01_78350_78851:Glyma.01G000322
#11      Gm01_94123_94624        Glyma.01G000359 0       gACR1   cds     Gm01_94123_94624:Glyma.01G000359
#12      Gm01_99690_100191       Glyma.01G000400 0       gACR1   5utr    Gm01_99690_100191:Glyma.01G000400
#13      Gm01_111474_111975      Glyma.01G000600 -8444   dACR    dACR    Gm01_111474_111975:Glyma.01G000600



def step01_build_acr_pairs_fl (input_acr_fl,input_gene_gff_fl,input_output_dir):

    step01_build_acr_pairs_fl_dir = input_output_dir + '/step01_build_acr_pairs_fl_dir'
    if not os.path.exists(step01_build_acr_pairs_fl_dir):
        os.makedirs(step01_build_acr_pairs_fl_dir)


    store_final_line_list = []
    with open (input_gene_gff_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if not eachline.startswith('#'):

                if col[2] == 'gene':

                    #if re.match('ID=.+',col[8]):
                    mt = re.match('ID=(.+)',col[8])
                    geneID = mt.group(1)

                    final_line = col[0] + '\t' + col[3] + '\t' + col[4] + '\t' + geneID + '\t' + col[5] + '\t' + col[6]
                    store_final_line_list.append(final_line)

    with open (step01_build_acr_pairs_fl_dir + '/temp_gene.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + step01_build_acr_pairs_fl_dir + '/temp_gene.bed > ' + step01_build_acr_pairs_fl_dir + '/temp_gene.sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    with open (input_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if not col[0].startswith('chloro'):
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[0] + '_' + col[1] + '_' + col[2]
                store_final_line_list.append(final_line)

    with open (step01_build_acr_pairs_fl_dir + '/temp_acr.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + step01_build_acr_pairs_fl_dir + '/temp_acr.bed > ' + step01_build_acr_pairs_fl_dir + '/temp_acr.sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools closest -a ' + step01_build_acr_pairs_fl_dir + '/temp_acr.sorted.bed' + \
          ' -b ' + step01_build_acr_pairs_fl_dir + '/temp_gene.sorted.bed' + \
          ' -d > ' + step01_build_acr_pairs_fl_dir + '/temp_acr_gene_closest.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_final_line_list = []
    first_line = 'peak_id' + '\t' +  'gene_id.close' + '\t' + 'dis2gene' + '\t' + 'peak.type1' + '\t' + 'peak.type2' + '\t' + 'peak2gene_id'
    store_final_line_list.append(first_line)
    count = 0
    with open (step01_build_acr_pairs_fl_dir + '/temp_acr_gene_closest.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col =  eachline.strip().split()
            count += 1

            final_line = str(count) + '\t' + col[3] + '\t' + col[7] + '\t' + col[10] + '\t' + 'na' + '\t' + 'na' + '\t' + col[3] + ':' + col[7]
            store_final_line_list.append(final_line)

    with open (step01_build_acr_pairs_fl_dir + '/opt.ACR.sorted_acr2gene.annotated.uniqACR.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def step02_getPerM_scRNAseq (input_required_script_dir,input_output_dir,input_scrna_object_fl):

    step02_getPerM_scRNAseq_dir = input_output_dir + '/step02_getPerM_scRNAseq_dir'
    if not os.path.exists(step02_getPerM_scRNAseq_dir):
        os.makedirs(step02_getPerM_scRNAseq_dir)

    ipt_target_colnm = 'new_cell_type'

    cmd = 'Rscript ' + input_required_script_dir + '/getPerM_celltypes.R' + \
          ' ' + input_scrna_object_fl + \
          ' ' + ipt_target_colnm + \
          ' ' + step02_getPerM_scRNAseq_dir
    print(cmd)
    subprocess.call(cmd,shell=True)


def step03_run_peak2gene (input_required_script_dir,input_acr_cpm_fl,input_core_num,input_output_dir):

    step03_run_peak2gene_dir = input_output_dir + '/step03_run_peak2gene_dir'
    if not os.path.exists(step03_run_peak2gene_dir):
        os.makedirs(step03_run_peak2gene_dir)

    ipt_target_colnm = 'new_cell_type'

    input_gene_cpm_fl =  input_output_dir + '/step02_getPerM_scRNAseq_dir/' + '/opt_perM_peaks_accessibility_' + ipt_target_colnm + ".txt"

    ##the first step is to match the pmat and gmat
    cmd = 'Rscript ' + input_required_script_dir + '/match_two_datasets.R' + \
          ' ' + input_acr_cpm_fl + \
          ' ' + input_gene_cpm_fl + \
          ' ' + step03_run_peak2gene_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

    ipt_acr_pairs_fl =  input_output_dir + '/step01_build_acr_pairs_fl_dir/temp_acr_gene_closest.txt'
    pmat = step03_run_peak2gene_dir + '/opt_pmat.txt'
    gmat = step03_run_peak2gene_dir + '/opt_gmat.txt'
    prefix = "pearl_millet_celltype"
    cmd = 'Rscript ' + input_required_script_dir + '/01_ACR2gene_atlas_celltype.R' + \
          ' ' + ipt_acr_pairs_fl + \
          ' ' + input_core_num + \
          ' ' + pmat + \
          ' ' + gmat + \
          ' ' + prefix + \
          ' ' + step03_run_peak2gene_dir + \
          ' ' + input_required_script_dir
    print(cmd)
    subprocess.call(cmd,shell=True)


    ##link the genes
    ipt_acr2gene = input_output_dir + '/' + 'step01_build_acr_pairs_fl_dir' + '/opt.ACR.sorted_acr2gene.annotated.uniqACR.bed'

    ipt_fdr_cutoff = '1'
    ipt_pval_cutoff = '0.05'

    cmd = 'Rscript ' + input_required_script_dir + '/02_Check_acr2gene_links_v5.R' + \
          ' ' + step03_run_peak2gene_dir + '/' + prefix + '_ACR2gene_pairs.observed.gene-peak-links.txt' + \
          ' ' + ipt_acr2gene + \
          ' ' + prefix + \
          ' ' + step03_run_peak2gene_dir + \
          ' ' + input_required_script_dir + \
          ' ' + ipt_fdr_cutoff + \
          ' ' + ipt_pval_cutoff
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##add the number of ACRs for the nearby genes and distance
    store_gene_acrinfor_dic = {}
    count = 0
    with open (ipt_acr2gene,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                geneID = col[2]
                peakID = col[1]
                dis2gene = col[3]
                peakIDdis2gene = peakID + ':' + dis2gene

                if geneID in store_gene_acrinfor_dic:
                    store_gene_acrinfor_dic[geneID].append(peakIDdis2gene)
                else:
                    store_gene_acrinfor_dic[geneID] = []
                    store_gene_acrinfor_dic[geneID].append(peakIDdis2gene)

    ipt_gene_link_annot_fl = input_output_dir + '/step03_run_peak2gene_dir/' + prefix + '.observed.gene-peak-links.annotated.txt'
    store_final_line_list = []
    count = 0
    with open (ipt_gene_link_annot_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                geneID = col[2]
                peakInfor = store_gene_acrinfor_dic[geneID]

                peakInfor_str = ';'.join(peakInfor)
                peakNum = len(peakInfor)

                final_line = eachline + '\t' + peakInfor_str + '\t' + str(peakNum)
                store_final_line_list.append(final_line)

            else:

                final_line = eachline + '\t' + 'all_peak_dis2gene' + '\t' + 'surround_peak_num'
                store_final_line_list.append(final_line)

    with open (step03_run_peak2gene_dir + '/' + prefix + '.observed.gene-peak-links.annotated.allpeakInfor.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')








def step04_select_top_gene_to_plot (input_output_dir,input_atac_meta_fl,input_select_pair_id_list_fl,input_required_script_dir):

    step04_select_top_gene_to_plot_dir = input_output_dir + '/step04_select_top_gene_to_plot_dir'
    if not os.path.exists(step04_select_top_gene_to_plot_dir):
        os.makedirs(step04_select_top_gene_to_plot_dir)

    pmat = input_output_dir + '/step03_run_peak2gene_dir' + '/opt_pmat.txt'
    gmat = input_output_dir + '/step03_run_peak2gene_dir' + '/opt_gmat.txt'

    prefix = "pearl_millet_celltype"
    ipt_gene_link_annot_fl = input_output_dir + '/step03_run_peak2gene_dir/' + prefix + '.observed.gene-peak-links.annotated.txt'

    ipt_R_script = input_required_script_dir + '/plot_corr_target_links.R'
    ipt_top_link_num = '10'
    ipt_type_str = 'pLink'


    cmd = 'Rscript ' + ipt_R_script + \
          ' ' + pmat + \
          ' ' + gmat + \
          ' ' + ipt_gene_link_annot_fl + \
          ' ' + ipt_top_link_num + \
          ' ' + step04_select_top_gene_to_plot_dir + \
          ' ' + ipt_type_str + \
          ' ' + input_atac_meta_fl + \
          ' ' + input_select_pair_id_list_fl
    print(cmd)
    subprocess.call(cmd,shell=True)






step01_build_acr_pairs_fl (input_acr_fl,input_gene_gff_fl,input_output_dir)
step02_getPerM_scRNAseq (input_required_script_dir,input_output_dir,input_scrna_object_fl)
step03_run_peak2gene (input_required_script_dir,input_acr_cpm_fl,input_core_num,input_output_dir)
#step04_select_top_gene_to_plot (input_output_dir,input_required_script_dir)




























