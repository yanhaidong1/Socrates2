#!/usr/bin/env python

import argparse
import glob
import sys
import os
import subprocess

##this pipeline is to load the data and perform the QC

def get_parsed_args():

    parser = argparse.ArgumentParser(description="Load data and perform QC")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    ##Users could choose to run all the QC at one time
    ##Or finish the step01 to obain the object and perform the following findcells or is cells by trying different parameters

    ##################################
    ##step01 obtain the Socrate object
    parser.add_argument("-open_build_obj", dest = 'open_build_object', help = 'Build an object first to store all information.'
                                                                              'Default: yes')

    parser.add_argument("-tn5_fl", dest = 'tn5_file', help='Users will provide a BED file with the Tn5 integration sites.')

    parser.add_argument("-Ggff_fl", dest='gene_gff_file',
                        help='Provide the gene gff file.')

    parser.add_argument("-chrsize_fl", dest = 'chrom_size_file', help='Provide a chrom size file.')


    ##parameter settings (required)
    parser.add_argument("-Pt_chr_name", dest = 'Pt_chrom_name',help = 'Provide a Pt chrom name within genome.')

    parser.add_argument("-Mt_chr_name", dest = 'Mt_chrom_name', help = 'Provide a Mt chrom name within genome.')


    ##parameter settings (optional)
    parser.add_argument("-prefix", dest='prefix_name', help='Provide a prefix added in the output files.'
                                                            'Default: output')

    parser.add_argument("-shift", dest = 'macs2_shift', help = 'Provide a parameter settings in MACS2.'
                                                               'Default: -75')

    parser.add_argument("-extsize", dest = 'macs2_extsize', help = 'Provide a parameter setting in MACS2.'
                                                                   'Default: 150')

    parser.add_argument("-fdr", dest = 'macs2_fdr', help = 'Provide a parameter setting to filter peak in MACS2.'
                                                           'Default: 0.1')

    parser.add_argument("-tss_size", dest = 'tss_window_size', help = 'A promoter TSS window size defined to perform the QC.'
                                                                      'Default: 2000')

    ###################
    ##step02 find cells
    parser.add_argument("-open_find_cells", dest = 'open_find_cells_step', help = 'Perform intital QC to find an overall data quality.'
                                                                                  'Default: yes')

    parser.add_argument("-ipt_raw_obj", dest = 'input_raw_object', help = 'Provide an object file obtained from the -open_build_obj.')


    parser.add_argument("-min_cells", dest = 'min_cell_val', help = 'Lower limit on the number of identified cells.'
                                                                    'Defaults: 1000')

    parser.add_argument("-max_cells", dest = 'max_cell_val', help = 'Upper limit on the number of identified cells.'
                                                                    'Defaults: 15000')

    parser.add_argument("-min_tn5", dest = 'min_tn5_val', help = 'Lower threshold for the minimum number of Tn5 integration sites for retaining obj barcode.'
                                                                 'Default: 1000')

    parser.add_argument("-org_flt_thresh", dest = 'organelle_filter_cutoff', help = 'Remove cells with an organelle ratio (Organalle/Total_reads) greater than N.'
                                                                                    'Default: 0.8')

    parser.add_argument("-tss_min_freq", dest = 'tss_min_freq_val', help = 'Minimum frequency of Tn5 sites near TSSs.'
                                                                           'Default: 0.2')

    parser.add_argument("-tss_z_thresh", dest = 'tss_z_thresh_val', help = 'Z-score threshold to remove barcodes below the mean.'
                                                                           'Default: 2')

    parser.add_argument("-frip_min_freq", dest = 'frip_min_freq_val', help = 'Minimum frequency of Tn5 sites within ACRs.'
                                                                           'Default: 0.2')

    parser.add_argument("-frip_z_thresh", dest = 'frip_z_thresh_val', help = 'Z-score threshold to remove barcodes with values below X standard deviations from the mean.'
                                                                             'Default: 2')

    ################
    ##step03 iscells
    parser.add_argument("-open_is_cells", dest = 'open_is_cells_step', help = 'This function compares each barcode to obj background and cell bulk reference to identify barcodes representing ambient DNA or broken nuclei.'
                                                                              'Default: yes')

    parser.add_argument("-ipt_findcell_obj", dest = 'input_findcell_object', help = 'Provide an object obtained from the -open_find_cells.')

    parser.add_argument("-num_test", dest='num_test_val',
                        help='Number of barcodes to query (ranked by total # of unique Tn5 insertions).'
                             'Default: 20000')

    parser.add_argument("-num_tn5", dest = 'num_tn5_val',help = 'Set the minimum number of Tn5 insertions to select test cells. Overridden by num.test.'
                                                                'Default: NULL')


    parser.add_argument("-num_REF", dest = 'num_ref_val',help = 'Number of cells to use as the cell bulk reference (top X cells based on # unique Tn5 insertions)'
                                                                'Default: 1000')

    parser.add_argument("-background_cutoff", dest = 'background_cutoff_val',help = 'Maximum unique Tn5 insertions to use for selecting barcodes for the background reference set)'
                                                                'Default: 100')


    #######################
    ##generate final matrix
    parser.add_argument("-win_size", dest = 'window_size', help = 'Divide the genome into different bins/windows of each being a specific size/bp. '
                                                                  'Default: 500')


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

    input_required_scripts_dir = args.required_script_dir

    store_final_parameter_line_list = []

    ##set the open of big steps
    if args.open_build_object is None:
        open_build_object_final = 'yes'
    else:
        if args.open_build_object == 'yes':
            open_build_object_final = 'yes'
        else:
            open_build_object_final = 'no'

    if args.open_find_cells_step is None:
        open_find_cells_step_final = 'yes'
    else:
        if args.open_find_cells_step == 'yes':
            open_find_cells_step_final = 'yes'
        else:
            open_find_cells_step_final = 'no'

    if args.open_is_cells_step is None:
        open_is_cells_step_final = 'yes'
    else:
        if args.open_is_cells_step == 'yes':
            open_is_cells_step_final = 'yes'
        else:
            open_is_cells_step_final = 'no'

    ##set the required files
    if open_build_object_final == 'yes':

        if args.tn5_file is None:
            print('Cannot find tn5 file, please provide it')
            return
        else:
            try:
                file = open(args.tn5_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the tn5 file!')
                return

        if args.gene_gff_file is None:
            print('Cannot find gene gff file, please provide it')
            return
        else:
            try:
                file = open(args.gene_gff_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the gene gff file!')
                return

        if args.chrom_size_file is None:
            print('Cannot find chromosome size file, please provide it')
            return
        else:
            try:
                file = open(args.chrom_size_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the chromosome size file!')
                return


    ##only if the build object do not open and find cells open we will provide the input object
    if open_build_object_final == 'no':

        if open_find_cells_step_final == 'yes':

            if args.input_raw_object is None:
                print('Cannot find a raw built object file, please provide it')
                return
            else:
                try:
                    file = open(args.input_raw_object, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the raw built object file!')
                    return

        else:

            if open_is_cells_step_final == 'yes':

                if args.input_findcell_object is None:
                    print('Cannot find a find cell object file, please provide it')
                    return
                else:
                    try:
                        file = open(args.input_findcell_object, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the find cell object file!')
                        return


    ##required parameters
    if args.Pt_chrom_name is not None:
        Pt_chrom_name_final = args.Pt_chrom_name
    else:
        print('Please provide the pt chromosome name')
        return

    if args.Mt_chrom_name is not None:
        Mt_chrom_name_final = args.Mt_chrom_name
    else:
        print('Please provide the mt chromosome name')
        return

    store_final_parameter_line_list.append('Pt <- ' + args.Pt_chrom_name)
    store_final_parameter_line_list.append('Mt <- ' + args.Mt_chrom_name)

    ##step01 parameters
    if args.prefix_name is not None:
        prefix = args.prefix_name
    else:
        prefix = 'output'

    store_final_parameter_line_list.append('out <- ' + prefix)


    if args.macs2_shift is not None:
        macs2_shift_final = args.macs2_shift
    else:
        macs2_shift_final = '-75'

    store_final_parameter_line_list.append('macs2_shift_final <- ' + macs2_shift_final)

    if args.macs2_extsize is not None:
        macs2_extsize_final = args.macs2_extsize
    else:
        macs2_extsize_final = '150'

    store_final_parameter_line_list.append('macs2_extsize_final <- ' + macs2_extsize_final)

    if args.macs2_fdr is not None:
        macs2_fdr_final = args.macs2_fdr
    else:
        macs2_fdr_final = '0.1'

    store_final_parameter_line_list.append('macs2_fdr_final <- ' + macs2_fdr_final)


    if args.tss_window_size is not None:
        tss_window_size_final = args.tss_window_size
    else:
        tss_window_size_final = '2000'

    store_final_parameter_line_list.append('tss_window_size_final <- ' + tss_window_size_final)

    ##step02 parameters
    if args.min_cell_val is not None:
        min_cell_val_final = args.min_cell_val
    else:
        min_cell_val_final = '1000'

    store_final_parameter_line_list.append('min_cell_val_final <- ' + min_cell_val_final)

    if args.max_cells is not None:
        max_cell_val_final = args.max_cells
    else:
        max_cell_val_final = '16000'

    store_final_parameter_line_list.append('max_cell_val_final <- ' + max_cell_val_final)

    if args.min_tn5 is not None:
        min_tn5_val = args.min_tn5
    else:
        min_tn5_val = '1000'

    store_final_parameter_line_list.append('min_tn5_val <- ' + min_tn5_val)

    if args.org_flt_thresh is not None:
        organelle_filter_cutoff_final = args.org_flt_thresh
    else:
        organelle_filter_cutoff_final = '0.8'

    store_final_parameter_line_list.append('organelle_filter_cutoff_final <- ' + organelle_filter_cutoff_final)

    if args.tss_min_freq is not None:
        tss_min_freq_val_final = args.tss_min_freq
    else:
        tss_min_freq_val_final = '0.2'

    store_final_parameter_line_list.append('tss_min_freq_val_final <- ' + tss_min_freq_val_final)

    if args.tss_z_thresh is not None:
        tss_z_thresh_final = args.tss_z_thresh
    else:
        tss_z_thresh_final = '2'

    store_final_parameter_line_list.append('tss_z_thresh_final <- ' + tss_z_thresh_final)

    if args.frip_min_freq is not None:
        frip_min_freq_final = args.frip_min_freq
    else:
        frip_min_freq_final = '0.2'

    store_final_parameter_line_list.append('frip_min_freq_final <- ' + frip_min_freq_final)

    if args.frip_z_thresh is not None:
        frip_z_thresh_final = args.frip_z_thresh
    else:
        frip_z_thresh_final = '1'

    store_final_parameter_line_list.append('frip_z_thresh_final <- ' + frip_z_thresh_final)

    ##step03
    if args.num_test_val is not None:
        num_test_val_final = args.num_test_val
    else:
        num_test_val_final = '20000'

    store_final_parameter_line_list.append('num_test_val_final <- ' + num_test_val_final)

    if args.num_tn5_val is not None:
        num_tn5_val_final = args.num_tn5_val
    else:
        num_tn5_val_final = 'NULL'

    store_final_parameter_line_list.append('num_tn5_val_final <- ' + num_tn5_val_final)

    if args.num_ref_val is not None:
        num_ref_val_final = args.num_ref_val
    else:
        num_ref_val_final = '1000'

    store_final_parameter_line_list.append('num_ref_val_final <- ' + num_ref_val_final)

    if args.background_cutoff_val is not None:
        background_cutoff_val_final = args.background_cutoff_val
    else:
        background_cutoff_val_final = '100'

    store_final_parameter_line_list.append('background_cutoff_val_final <- ' + background_cutoff_val_final)

    if args.window_size is not None:
        window_size_bp = args.window_size
    else:
        window_size_bp = '500'

    store_final_parameter_line_list.append('win_size <- ' + window_size_bp)


    chr_size = 0
    with open (args.chrom_size_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            chr_size += int(col[1])

    store_final_parameter_line_list.append('genomesize <- ' + str(chr_size))

    ##we will firstly build up the parameter setting file
    with open (output_dir + '/temp_defined_parameters.config','w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')

    ipt_script = input_required_scripts_dir + '/utils_load_and_QC_data/load_and_QC.R'
    ipt_tn5_fl = args.tn5_file
    ipt_chr_size_fl = args.chrom_size_file
    ipt_gene_gff_fl = args.gene_gff_file
    ipt_config_fl = output_dir + '/temp_defined_parameters.config'
    ipt_path_to_preload_R = input_required_scripts_dir + '/preload_R'

    ##running
    cmd = 'Rscript ' + ipt_script + \
          ' ' + ipt_tn5_fl + \
          ' ' + ipt_chr_size_fl + \
          ' ' + ipt_gene_gff_fl + \
          ' ' + ipt_config_fl + \
          ' ' + output_dir + \
          ' ' + open_build_object_final + \
          ' ' + open_find_cells_step_final + \
          ' ' + open_is_cells_step_final + \
          ' ' + ipt_path_to_preload_R
    print(cmd)
    subprocess.call(cmd,shell=True)



if __name__ == "__main__":
    main()

















