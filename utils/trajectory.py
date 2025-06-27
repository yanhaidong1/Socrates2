#!/usr/bin/env python

import argparse
import glob
import sys
import subprocess
import os
import re



def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')


    parser.add_argument("-feature" ,dest = 'category_run', help= 'Users must provide a category to identify significant Gene/TF/Motif/ACR during cell progression.'
                                                                  'Default: Gene,TF,Motif,ACR.')

    #parser.add_argument("-ct_colnm", dest='celltype_cluster_col_name',
    #                    help='Define the cluster column name of the meta file.')


    ##############
    ##for the gene, it must have gene accessibility file
    parser.add_argument("-soc_obj_gene_acc", dest='soc_obj_gene_acc_fl',
                        help='Provide an object obtained from the Calculate gene accessibility step.')

    ##############
    ##for the peak, it must have the peak sparse file
    parser.add_argument("-ACR_sparse_fl",dest='peak_sparse_file',help='Provide a peak sparse file with Tn5 integration site.')

    ###############
    ##for the motif, it must provide the motif deviation score file
    parser.add_argument("-motif_dev_fl", dest = 'motif_deviation_file', help='Provide the motif deviation file with smoothed version.')

    ############
    ##for the TF
    parser.add_argument("-TF_fl", dest = 'TF_list_file',help = 'Provide a list of TF file')


    #######################
    ##provide the meta file
    parser.add_argument("-meta_fl", dest = 'meta_file', help = 'Provide a meta file that show the final cell identity. make sure the cell identity column name is cell_identity and celltype_color')

    ######################
    ##provide the svd file soc.obj contains the svd
    parser.add_argument("-soc_obj_svd", dest = 'soc_obj_svd_fl', help = 'Provide a svd file.')

    #############################
    ##provide the trajectory path
    parser.add_argument("-celltype_path", dest='celltype_path_str', help = 'Provide a cell type path string.')


    ##Optional parameters
    parser.add_argument("-core", dest='core_number', help='Specify how many cores we will use.'
                                                          'Default: 1')

    parser.add_argument('-prefix', dest = 'prefix_name', help = 'Specify a prefix name for the final output file.'
                                                                'Default: output')

    #parser.add_argument("-color_colnm", dest='color_cluster_col_name',
    #                    help='Define the color cluster column name of the meta file.')

    parser.add_argument("-TopACR_num",dest = 'Top_ACR_number_plot' ,help = 'number of top ACR will be plotted during trajectory analysis.'
                                              'Default: 10000')

    parser.add_argument("-open_equal_cellnum",dest = 'open_equal_cell_number', help = 'if users would like to initiate similar cell number per cell type.'
                                                                                      'Default:no')


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

    #if args.celltype_cluster_col_name is None:
    #    print('Please provide the name of column specifying the cell type identity in the meta file.')
    #    return

    input_required_scripts_dir = args.required_script_dir

    ##check if the category works
    if args.category_run is None:

        print('Users do not define the categories to run, and the default is to perform trajectory analysis for Gene, TF, Motif, and ACR.')

        final_run_category = 'Gene,TF,Motif,ACR'

    else:

        final_run_category = args.category_run

    if args.meta_file is None:
        print('Cannot find the meat file, please provide it')
        return
    else:
        try:
            file = open(args.meta_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the meta file!')
            return

    if args.soc_obj_svd_fl is None:
        print('Cannot find a object storing the svd file, please provide it')
        return
    else:
        try:
            file = open(args.soc_obj_svd_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the object storing the svd file!')
            return


    ##we will build a file to store all the information as the configure file
    store_final_parameter_line_list = []

    if args.celltype_path_str is None:
        print ('Users need to provide the cell type path. For example: -celltype_path celltype1,celltype2,celltype3')
        return
    else:
        celltype_path_final = args.celltype_path_str

        celltype_path_final_list = celltype_path_final.split(',')

        store_new_path_list = []
        for eachitem in celltype_path_final_list:
            final_item = '\'' + eachitem + '\''
            store_new_path_list.append(final_item)

        store_new_path_str = ','.join(store_new_path_list)
        traj = 'traj <- c(' + store_new_path_str + ')'
        store_final_parameter_line_list.append(traj)
        LC =  'LC <- c(' + store_new_path_str + ')'
        store_final_parameter_line_list.append(LC)


    ##after specifying the category, we need to check which category we will focuse on
    final_run_category_list = final_run_category.split(',')

    if 'Gene' in final_run_category_list:

        store_final_parameter_line_list.append('openGN <- ' + '\'' + 'yes' + '\'')

        if args.soc_obj_gene_acc_fl is None:
            print('Cannot find an object storing the gene accessibility, please provide it')
            return
        else:
            try:
                file = open(args.soc_obj_gene_acc_fl, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the object storing the gene accessibility file!')
                return
    else:
        store_final_parameter_line_list.append('openGN <- ' + '\'' + 'no' + '\'')

    if 'TF' in final_run_category_list:

        store_final_parameter_line_list.append('openTF <- ' + '\'' + 'yes' + '\'')

        if args.TF_list_file is None:
            print('Cannot find a TF list file, please provide it')
            return
        else:
            try:
                file = open(args.TF_list_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the TF list file!')
                return

        if args.soc_obj_gene_acc_fl is None:
            print('Cannot find an object storing the gene accessibility, please provide it')
            return
        else:
            try:
                file = open(args.soc_obj_gene_acc_fl, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the object storing the gene accessibility file!')
                return
    else:
        store_final_parameter_line_list.append('openTF <- ' + '\'' + 'no' + '\'')

    if 'Motif' in final_run_category_list:

        store_final_parameter_line_list.append('openMT <- ' + '\'' + 'yes' + '\'')

        if args.motif_deviation_file is None:
            print('Cannot find a motif deviation file, please provide it')
            return
        else:
            try:
                file = open(args.motif_deviation_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the motif deviation file!')
                return

    else:
        store_final_parameter_line_list.append('openMT <- ' + '\'' + 'no' + '\'')

    if 'ACR' in final_run_category_list:

        store_final_parameter_line_list.append('openACR <- ' + '\'' + 'yes' + '\'')

        if args.peak_sparse_file is None:
            print('Cannot find a peak sparse file, please provide it')
            return
        else:
            try:
                file = open(args.peak_sparse_file, 'r')  ##check if the file is not the right file
            except IOError:
                print('There was an error opening the peak sparse file!')
                return

    else:
        store_final_parameter_line_list.append('openACR <- ' + '\'' + 'no' + '\'')


    ##we will not open the TGN
    store_final_parameter_line_list.append('openTGN <- ' + '\'' + 'no' + '\'')

    with open(output_dir + '/temp_defined_parameters.config', 'w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')


    if args.prefix_name is not None:
        prefix_name_final = args.prefix_name
    else:
        prefix_name_final = 'output'


    if args.core_number is not None:
        core_number_run = args.core_number
    else:
        core_number_run = '1'


    #if args.color_cluster_col_name is not None:
    #    color_cluster_col_name_final = args.color_cluster_col_name
    #else:
    #    color_cluster_col_name_final = 'na'


    if args.Top_ACR_number_plot is not None:
        Top_ACR_number_plot_final = args.Top_ACR_number_plot
    else:
        Top_ACR_number_plot_final = '10000'

    if args.open_equal_cell_number is None:
        open_equal_cell_number_final = 'no'
        print('Users do not want to allow the cell type with equal number of cell')
    else:
        if args.open_equal_cell_number == 'yes':
            open_equal_cell_number_final = 'yes'
        else:
            open_equal_cell_number_final = 'no'
            print('Users do not want to allow the cell type with equal number of cell')



    ##############################################
    ##prepare the data according to the categories
    step01_prepare_sparse_fl_dir = output_dir + '/step01_prepare_sparse_fl_dir'
    if not os.path.exists(step01_prepare_sparse_fl_dir):
        os.makedirs(step01_prepare_sparse_fl_dir)

    ##if TF is in the list
    if 'TF' in final_run_category_list:

        print('- open prepare TF sparse')

        input_TF_fl = args.TF_list_file

        ipt_script = input_required_scripts_dir + '/utils_trajectory/build_TF_acc.R'
        ipt_gene_acc_rds_fl =  args.soc_obj_gene_acc_fl

        cmd = 'Rscript ' + ipt_script + \
              ' ' + ipt_gene_acc_rds_fl + \
              ' ' + input_TF_fl + \
              ' ' + step01_prepare_sparse_fl_dir
        print(cmd)
        subprocess.call(cmd,shell=True)

    if 'ACR' in final_run_category_list:

        print('- open prepare peak mtx')

        ipt_script = input_required_scripts_dir + '/utils_trajectory/transfer_to_mtx.R'
        ipt_peak_sparse_fl = args.peak_sparse_file

        cmd = 'Rscript ' + ipt_script + \
              ' ' + ipt_peak_sparse_fl + \
              ' ' + step01_prepare_sparse_fl_dir + \
              ' ' + 'peak_sparse.mtx'
        print(cmd)
        subprocess.call(cmd, shell=True)

    step02_running_trajectory_analysis_dir = output_dir + '/step02_running_trajectory_analysis_dir'
    if not os.path.exists(step02_running_trajectory_analysis_dir):
        os.makedirs(step02_running_trajectory_analysis_dir)

    ipt_script = input_required_scripts_dir + '/utils_trajectory/pseudotime.R'

    target_config_fl = output_dir + '/temp_defined_parameters.config'

    if 'TF' in final_run_category_list:
        target_tf_acc_threeCol_fl = step01_prepare_sparse_fl_dir + '/TF.acc.rds'
    else:
        target_tf_acc_threeCol_fl = 'na'

    if 'Motif' in final_run_category_list:
        target_motif_fl = args.motif_deviation_file
    else:
        target_motif_fl = 'na'


    target_meta_fl = args.meta_file
    target_svd_obj_fl = args.soc_obj_svd_fl

    if 'ACR' in final_run_category_list:
        target_peak_mtx_rds_fl = step01_prepare_sparse_fl_dir + '/' + 'peak_sparse.mtx' + '.rds'
    else:
        target_peak_mtx_rds_fl = 'na'

    if 'Gene' in final_run_category_list:
        target_allgene_acc_threCol_fl = args.soc_obj_gene_acc_fl
    else:
        target_allgene_acc_threCol_fl = 'na'



    ##step02_s1_sub_open_provided_reduced_meta_fl we will allow it to be the file path
    ##or 'na'
    step02_s1_sub_open_provided_reduced_meta_fl = 'na'
    input_target_gene_list_fl = 'na'


    ##we will design scripts to plot with the color
    ##Or we will add the color column besides the option 5 generate the final UMAP plot with the color column
    ##Or it will automatically generate several colors

    step02_s1_sub_target_cluster_col = 'cell_identity'
    step02_s1_sub_target_cluster_color = 'celltype_color'
    step02_s1_sub_plot_top_ACR_num = Top_ACR_number_plot_final
    step02_s1_sub_open_equal_sample_ctnum = open_equal_cell_number_final


    cmd = 'Rscript ' + ipt_script + \
          ' ' + target_peak_mtx_rds_fl + \
          ' ' + target_motif_fl + \
          ' ' + target_allgene_acc_threCol_fl + \
          ' ' + target_tf_acc_threeCol_fl + \
          ' ' + target_meta_fl + \
          ' ' + target_svd_obj_fl + \
          ' ' + prefix_name_final + \
          ' ' + target_config_fl + \
          ' ' + core_number_run + \
          ' ' + step02_running_trajectory_analysis_dir + \
          ' ' + step02_s1_sub_target_cluster_col + \
          ' ' + step02_s1_sub_target_cluster_color + \
          ' ' + input_target_gene_list_fl + \
          ' ' + step02_s1_sub_plot_top_ACR_num + \
          ' ' + step02_s1_sub_open_equal_sample_ctnum + \
          ' ' + step02_s1_sub_open_provided_reduced_meta_fl
    print(cmd)
    subprocess.call(cmd, shell=True)







if __name__ == "__main__":
    main()








































































