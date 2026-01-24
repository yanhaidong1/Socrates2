#!/usr/bin/env python

import argparse
#import glob
import sys
import subprocess
#import os
#import re


##this pipeline is to predict the activating and repressing ACRs

def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')


    ##########
    ##required directories for all the steps
    parser.add_argument("-soc_obj", dest='soc_object_fl',
                        help='Provide an object obtained from the call peak per cell type.')

    parser.add_argument("-cpm_fl", dest='peak_cpm_file',
                        help='Users will provide the CPM for the peak per cell type file.')
    ##/scratch/hy17471/database_031724/section3_add_pearlmillet_101424/08_diff_peak_071325/output_dir_18kcells_cutoff0.01_updateAnnot/s1_open_prepare_peak_tn5_final/opt_peaks_celltype_CPM_dir/opt_perM_peaks_accessibility_cell_identity.txt

    parser.add_argument("-gene_gff", dest='gene_gff_file', help='Users must provide gene gff file.')

    parser.add_argument("-scrna_obj", dest='scrna_object_file',
                        help='Users must provide the scrna-seq object from obtained from Seurat tool.')


    ##########
    ##optional if users provide the acr_fl or tn5_fl seperately
    parser.add_argument("-acr_fl", dest = 'peak_file',help = 'Users must provide the acr file for the target species.')

    #parser.add_argument("-meta_fl", dest = 'meta_file', help = 'Users will provide the meta file.')

    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')




    ##parse of parameters
    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    if args.required_script_dir is None:
        print('Cannot find required script dir, please provide the dir in \'-script_dir\' !')
        return
    else:
        input_required_scripts_dir = args.required_script_dir


    if args.soc_object_fl is None:
        use_soc_object = 'no'
        print('Users choose to provide ACR and meta file')



        if args.peak_file is None:
            print('Cannot find the peak file, please provide it')
            return
        else:
            try:
                file = open(args.peak_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the peak file!')
                return

    else:
        use_soc_object = 'yes'

        try:
            file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the soc object file!')
            return


    if args.peak_cpm_file is None:
        print('Cannot find the cpm file, please provide it')
        return
    else:
        try:
            file = open(args.peak_cpm_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the cpm file!')
            return

    if args.gene_gff_file is None:
        print('Cannot find the gene gff file, please provide it')
        return
    else:
        try:
            file = open(args.gene_gff_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the gene gff file!')
            return

    if args.scrna_object_file is None:
        print('Cannot find the scrna object file, please provide it')
        return
    else:
        try:
            file = open(args.scrna_object_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the scrna object file!')
            return

    if args.core_number is not None:
        core_number_run = args.core_number
    else:
        core_number_run = '1'

    ##run the correlation
    if use_soc_object == 'yes':

        ipt_object_fl = args.soc_object_fl
        ##extract the final meta and peak file
        cmd = 'Rscript ' + input_required_scripts_dir + '/utils_predict_activating_repressing_ACR/read_object.R' + \
              ' ' + ipt_object_fl + \
              ' ' + output_dir
        print(cmd)
        subprocess.call(cmd,shell=True)

        ipt_peak_fl = output_dir + '/final_peak.txt'
        #ipt_meta_fl = output_dir + '/final_meta.txt'

    else:
        ipt_peak_fl = args.peak_file
        #ipt_meta_fl = args.meta_file

    ipt_gff_fl = args.gene_gff_file
    ipt_scrna_obj_fl = args.scrna_object_file
    ipt_cpm_fl = args.peak_cpm_file

    cmd = 'python ' + input_required_scripts_dir + '/utils_predict_activating_repressing_ACR/pipeline_correlation_ACR_gene.py' + \
          ' ' + input_required_scripts_dir + '/utils_predict_activating_repressing_ACR/input_required_scripts_dir' + \
          ' ' + ipt_peak_fl + \
          ' ' + ipt_gff_fl + \
          ' ' + ipt_scrna_obj_fl + \
          ' ' + ipt_cpm_fl + \
          ' ' + 'na' + \
          ' ' + core_number_run + \
          ' ' + output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)




if __name__ == "__main__":
    main()

