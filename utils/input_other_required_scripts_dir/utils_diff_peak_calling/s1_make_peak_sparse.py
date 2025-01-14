#!/usr/bin/env python

##updating 112222 we will build the data for the single organ
##updating 062122 we will use all the peaks other than the organ specific peaks to find the right
##updating 060222 set an option to call n binary peaks
##updating 052422 build the sparse with ACR for each chromosome
##updating 052322 we will make a combined version of sparse file
##updating 052022 we will generate the peak information
##updating 041522 we need to use the new peak file to generate the sparse file
##this script we will make the peak sparse file for each organ with using the parallele

##this script is to call peaks for each lib
import re
import glob
import subprocess
import os



#input_allorgans_peak_fl = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_filterPeak_040722/pipeline_ver/opt_final_collect_organ_peaks_combineAll_052322/opt_final_peak_sorted.txt

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_filterPeak_040722/pipeline_ver/opt_final_collect_organ_peaks_052022

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_filterPeak_040722/pipeline_ver/output_dir_multipart/store_all_final_organ_flt_peaks_dir_sixpart_threepartBase

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/output_dir_040322/opt_store_all_organs_unique_peak_dir/

#input_tn5_bed_fl = sys.argv[2]
#input_all_organ_tn5_dir = sys.argv[2]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/01_post_process_build_meta_pipeline_062521_011322/collect_all_tn5_bed_dir_choose_for_sep_analysis/output_dir_combine_sort_ver

#input_genome_fai_fl = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/Osativa323v7.fa.fai.sorted

#input_config_fl = sys.argv[3]
##

#input_fastSparsetn5_pl = sys.argv[4]
##fastSparse.tn5.pl

##updating 060222
##generate a non binary case
#input_nfastSparsetn5_pl = sys.argv[5]

#input_output_dir = sys.argv[5]




def make_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
                 input_output_dir):

    prefix = 'opt'

    opt_sorted_peak_bed_dir = input_output_dir + '/opt_sorted_peak_bed_dir'
    if not os.path.exists(opt_sorted_peak_bed_dir):
        os.makedirs(opt_sorted_peak_bed_dir)

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'
    if not os.path.exists(opt_peak_sparse_dir):
        os.makedirs(opt_peak_sparse_dir)

    peak_bed_fl_path = input_peak_fl


    #peak_bed_fl_path = eachdir + '/' + dir_nm + '.unique500bpPeaks.bed'
    sorted_tn5_bed_fl_path = input_tn5_bed_fl

    cmd = 'sort -k1,1V -k2,2n ' + peak_bed_fl_path + ' > ' + opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)
    sorted_peak_bed_fl_path = opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'

    cmd = 'bedtools intersect -a ' +  sorted_tn5_bed_fl_path + \
          ' -b ' + sorted_peak_bed_fl_path + ' -wa -wb' + \
          ' -sorted | perl ' + input_fastSparsetn5_pl + \
          ' - > ' + opt_peak_sparse_dir + '/opt_peak.sparse'
    print(cmd)
    subprocess.call(cmd,shell=True)


def sort_sparse (input_output_dir):

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'

    opt_peak_sparse_sorted_dir = input_output_dir + '/opt_peak_sparse_sorted_dir'
    if not os.path.exists(opt_peak_sparse_sorted_dir):
        os.makedirs(opt_peak_sparse_sorted_dir)

    allsparse_fl_list = glob.glob(opt_peak_sparse_dir + '/*')
    for eachfl in allsparse_fl_list:
        mt = re.match('.+/(.+)\.sparse',eachfl)
        flnm = mt.group(1)

        cmd = 'sort -k1,1V -k2,2n ' + eachfl + ' > ' + opt_peak_sparse_sorted_dir + '/' + flnm + '_sorted.sparse'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##remove the peak sparse file
        cmd = 'rm ' + eachfl
        print(cmd)
        subprocess.call(cmd,shell=True)







