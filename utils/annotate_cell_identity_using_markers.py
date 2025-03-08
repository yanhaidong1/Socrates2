#!/usr/bin/env python

import argparse
import glob
import sys
import os
import subprocess
import re

##updating 022025 add GO enrichment test
##this pipeline is to annotate cell identity including several methods:
##these methods are mainly based on the marker genes
##1) visualize marker genes using UMAP
##2) visualize marekr genes using dot plot
##3) aggregate the marker performance
##4) etc.

def get_parsed_args():

    parser = argparse.ArgumentParser(description="cell type peak calling pipeline")

    ##require files
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    ##required directories for all the steps
    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    ##step 01
    parser.add_argument("-open_markerUMAP", dest='open_vis_marker_UMAP',help = "Run visualization of marker genes using UMAP.")

    ##updating 010425
    parser.add_argument("-soc_obj", dest='soc_object_fl', help='Provide an object obtained from the Calculate gene accessibility step.')

    #parser.add_argument("-impute_gene_mtx", dest= 'impute_gene_acc_matrix', help = 'Provide the impute gene accessibility file.')

    #arser.add_argument("-meta_fl", dest='meta_file', help='Provide a meta file aftering the clustering.')

    parser.add_argument("-marker_fl", dest='marker_file',help = 'Provide a marker gene file.')

    ##Optional parameters
    parser.add_argument("-lim_val", dest = 'vis_lim_value', help = 'Provide a lim value of visualizing the markers in UMAP.')

    ##step 02
    parser.add_argument("-open_dotplot", dest = 'open_vis_marker_dotplot', help = 'Run visualization of marker genes using dotplot.')

    parser.add_argument("-gene_tn5_fl", dest='gene_tn5_sparse_file',
                        help='Provide file showing a gene with tn5 insertion per cell.')

    parser.add_argument("-clustnm", dest='target_cluster_name',
                        help='Provide a target cluster name within a cluster file.')

    parser.add_argument("-open_only_plot", dest = 'open_only_dotplot', help = 'If open only plot, users need to provide a rds file with data already prepared.'
                                                                              'Default: no')

    parser.add_argument("-updated_soc_obj", dest = 'updated_soc_object_fl', help = 'Once we open only plot, we need to specify a prepared gene file used for the dotplot.')


    parser.add_argument("-dotplotwidth", dest='dotplotwidth_val',
                        help='Specify the dot plot width value.'
                             'Default: 10')

    parser.add_argument("-dotplotheight", dest='dotplotheight_val',
                        help='Specify the dot plot height value.'
                             'Default: 10')

    ##step 03
    parser.add_argument("-open_aggregateAnnot", dest = 'open_aggregate_annotation', help = 'Run annotation based on aggregating marker genes.')

    parser.add_argument("-gene_acc_mtx", dest = 'gene_accessibility_mtx', help = 'Provide the gene accessibility matrix file.')

    parser.add_argument("-svd_fl", dest='svd_file', help='Provide a svd file aftering the clusteirng')

    ##step 04

    parser.add_argument("-open_GO_annot", dest = 'open_GO_enrich_annot', help = 'Run Go enrichment analysis for the group of genes enriched in a specific cell type.')

    parser.add_argument("-gene_GO_fl", dest = 'gene_GO_file', help = 'Provide a gene with GO process file that helps the GO enrichment test.')

    parser.add_argument("-log2fc_cutoff", dest = 'log2fc_cutoff_val', help = 'Set a threshold to filter the genes with log2(fold change) above a specific value.'
                                                                               'Default: 0.25.')




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



    if args.open_vis_marker_UMAP is None:

        open_markerUMAP_final = 'no'

    else:
        if args.open_vis_marker_UMAP == 'yes':

            if args.soc_object_fl is None:
                print('Cannot find the soc object file, please provide it')
                return
            else:
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the soc object file!')
                    return

            open_markerUMAP_final = 'yes'

            if args.marker_file is None:
                print('Cannot find marker gene file, please provide it')
                return
            else:
                try:
                    file = open(args.marker_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the marker gene file!')
                    return

        else:
            open_markerUMAP_final = 'no'
            print('Users choose to close the marker visualization in UMAP, please use \'-open_markerUMAP yes\' to open this step')


    open_only_dotplot_final = 'no'

    if args.open_vis_marker_dotplot is None:

        open_vis_marker_dotplot_final = 'no'

    else:
        if args.open_vis_marker_dotplot == 'yes':

            open_vis_marker_dotplot_final = 'yes'

            if args.open_only_dotplot is None:
                open_only_dotplot_final = 'no'
            else:
                if args.open_only_dotplot == 'yes':
                    open_only_dotplot_final = 'yes'
                else:
                    print('Please set -open_only_plot yes to open this step')
                    return


            if open_only_dotplot_final == 'no':

                if args.soc_object_fl is None:
                    print('Cannot find the soc object file, please provide it')
                    return
                else:
                    try:
                        file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the soc object file!')
                        return

                if args.marker_file is None:
                    print('Cannot find marker gene file, please provide it')
                    return
                else:
                    try:
                        file = open(args.marker_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the marker gene file!')
                        return

            else:

                if args.updated_soc_object_fl is None:
                    print('Cannot find updated soc object file, please provide it')
                    return
                else:
                    try:
                        file = open(args.updated_soc_object_fl, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening updated soc object file!')
                        return

                if args.marker_file is None:
                    print('Cannot find marker gene file, please provide it')
                    return
                else:
                    try:
                        file = open(args.marker_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the marker gene file!')
                        return


        else:
            open_vis_marker_dotplot_final = 'no'
            print(
                'Users choose to close the marker visualization based on dotplot, please use \'-open_dotplot yes\' to open this step')


    if args.open_aggregate_annotation is None:
        open_aggregate_annotation_final = 'no'
    else:

        if args.open_aggregate_annotation == 'yes':

            open_aggregate_annotation_final = 'yes'

            if args.soc_object_fl is None:
                print('Cannot find the soc object file, please provide it')
                return
            else:
                try:
                    file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the soc object file!')
                    return

            if args.marker_file is None:
                print('Cannot find marker gene file, please provide it')
                return
            else:
                try:
                    file = open(args.marker_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the marker gene file!')
                    return

        else:
            open_aggregate_annotation_final = 'no'



    ##updating 022025
    if args.open_GO_enrich_annot is None:
        open_GO_enrich_annot_final = 'no'
    else:
        if args.open_GO_enrich_annot == 'yes':

            open_GO_enrich_annot_final = 'yes'

            if args.gene_GO_file is None:
                print('Cannot find gene GO file, please provide it')
                return
            else:
                try:
                    file = open(args.gene_GO_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the gene GO file!')
                    return
        else:
            open_GO_enrich_annot_final = 'no'


    if args.log2fc_cutoff_val is not None:
        log2fc_cutoff_val_final = args.log2fc_cutoff_val
    else:
        log2fc_cutoff_val_final = '0.25'






    if args.vis_lim_value is not None:
        vis_lim_value_final = args.vis_lim_value
    else:
        vis_lim_value_final = '0.999'


    if args.target_cluster_name is not None:
        target_cluster_nm_final = args.target_cluster_name
    else:
        target_cluster_nm_final = 'LouvainClusters'

    if args.dotplotwidth_val is not None:
        dotplotwidth_val_final = args.dotplotwidth_val
    else:
        dotplotwidth_val_final = '10'

    if args.dotplotheight_val is not None:
        dotplotheight_val_final = args.dotplotheight_val
    else:
        dotplotheight_val_final = '10'




    if open_markerUMAP_final == 'yes':

        print('Users choose to visualize marker gene accessibility based on UMAP')

        ##updating 010425 we will get the input prefix of the soc obj fl
        input_soc_obj_fl = args.soc_object_fl
        if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
            mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
            input_prefix = mt.group(1)
        else:
            if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                input_prefix = mt.group(1)
            else:
                print('Please use *.atac.soc.rds file without changing the file name')
                return

        open_markerUMAP_final_dir = output_dir + '/open_markerUMAP_final_dir'
        if not os.path.exists(open_markerUMAP_final_dir):
            os.makedirs(open_markerUMAP_final_dir)


        vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_UMAP.R'

        ##updating 010524
        #ipt_meta_fl = args.meta_file
        #ipt_impute_fl = args.impute_gene_acc_matrix
        ipt_marker_fl = args.marker_file

        cmd = 'Rscript ' + vis_R_script + \
              ' ' + input_soc_obj_fl + \
              ' ' + ipt_marker_fl + \
              ' ' + vis_lim_value_final + \
              ' ' + open_markerUMAP_final_dir
        print(cmd)
        subprocess.call(cmd,shell=True)



    if open_vis_marker_dotplot_final == 'yes':

        print('Users choose to visualize marker gene accessibility based on dotplot')

        open_vis_marker_dotplot_final_dir = output_dir + '/open_vis_marker_dotplot_final_dir'
        if not os.path.exists(open_vis_marker_dotplot_final_dir):
            os.makedirs(open_vis_marker_dotplot_final_dir)

        if open_only_dotplot_final == 'no':

            vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_dotplot_prepare.R'

            ##updating 010425 we will get the input prefix of the soc obj fl
            input_soc_obj_fl = args.soc_object_fl
            if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                input_prefix = mt.group(1)
            else:
                if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                    mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                    input_prefix = mt.group(1)
                else:
                    print('Please use *.atac.soc.rds file without changing the file name')
                    return

            #ipt_impute_fl = args.impute_gene_acc_matrix
            #ipt_sparse_fl = args.gene_tn5_sparse_file
            #ipt_meta_fl = args.meta_file

            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + input_soc_obj_fl + \
                  ' ' + target_cluster_nm_final + \
                  ' ' + open_vis_marker_dotplot_final_dir + \
                  ' ' + input_prefix
            print(cmd)
            subprocess.call(cmd,shell=True)

            vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_dotplot_plotting.R'

            #ipt_prepared_gene_fl = output_dir + '/open_vis_marker_dotplot_final_dir/opt_gene_zscore_accprop.rds'
            ipt_marker_fl = args.marker_file

            input_soc_obj_fl = open_vis_marker_dotplot_final_dir + '/' + input_prefix + '.atac.soc.rds'

            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + input_soc_obj_fl + \
                  ' ' + ipt_marker_fl + \
                  ' ' + open_vis_marker_dotplot_final_dir + \
                  ' ' + dotplotwidth_val_final + \
                  ' ' + dotplotheight_val_final
            print(cmd)
            subprocess.call(cmd, shell=True)

        else:

            vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_dotplot_plotting.R'

            #ipt_prepared_gene_fl = output_dir + '/open_vis_marker_dotplot_final_dir/opt_gene_zscore_accprop.rds'
            ipt_marker_fl = args.marker_file

            ##updating 010425 we will get the input prefix of the soc obj fl
            input_soc_obj_fl = args.updated_soc_object_fl
            if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                input_prefix = mt.group(1)
            else:
                if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                    mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                    input_prefix = mt.group(1)
                else:
                    print('Please use *.atac.soc.rds file without changing the file name')
                    return


            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + input_soc_obj_fl + \
                  ' ' + ipt_marker_fl + \
                  ' ' + open_vis_marker_dotplot_final_dir + \
                  ' ' + dotplotwidth_val_final + \
                  ' ' + dotplotheight_val_final
            print(cmd)
            subprocess.call(cmd,shell=True)

    if open_aggregate_annotation_final == 'yes':

        print('Users choose to open the aggregating marker accessibility to visualize the markers')

        input_soc_obj_fl = args.soc_object_fl
        if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
            mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
            input_prefix = mt.group(1)
        else:
            if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                input_prefix = mt.group(1)
            else:
                print('Please use *.atac.soc.rds file without changing the file name')
                return

        open_aggregate_annotation_final_dir = output_dir + '/open_aggregate_annotation_final_dir'
        if not os.path.exists(open_aggregate_annotation_final_dir):
            os.makedirs(open_aggregate_annotation_final_dir)

        vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_aggregate_markers_UMAP.R'

        #ipt_gene_acc_mtx_rds_fl = args.gene_accessibility_mtx
        #ipt_meta_fl = args.meta_file
        #ipt_svd_fl = args.svd_file
        prefix = 'opt_annot'
        openAnnot = 'yes'
        ipt_marker_fl = args.marker_file

        cmd = 'Rscript ' + vis_R_script + \
              ' ' + input_soc_obj_fl + \
              ' ' + ipt_marker_fl + \
              ' ' + prefix + \
              ' ' + openAnnot + \
              ' ' + open_aggregate_annotation_final_dir + \
              ' ' + input_prefix
        print(cmd)
        subprocess.call(cmd,shell=True)


    ##udpating 022025
    if open_GO_enrich_annot_final == 'yes':

        print('Users choose to open the GO enrichment test for cell type genes')

        open_GO_enrich_annot_final_dir = output_dir + '/open_GO_enrich_annot_final_dir'
        if not os.path.exists(open_GO_enrich_annot_final_dir):
            os.makedirs(open_GO_enrich_annot_final_dir)

        ##step01 prepare the input of GO
        step01_prepare_GO_ipt_data_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/perform_GO_enrich_analysis.R'

        input_soc_obj_fl = args.soc_object_fl

        if re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl):
            mt = re.match('.+/(.+)\.atac\.soc\.rds', input_soc_obj_fl)
            input_prefix = mt.group(1)
        else:
            if re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl):
                mt = re.match('(.+)\.atac\.soc\.rds', input_soc_obj_fl)
                input_prefix = mt.group(1)
            else:
                print('Please use *.atac.soc.rds file without changing the file name')
                return

        input_gene_GO_fl = args.gene_GO_file

        cmd = 'python ' + step01_prepare_GO_ipt_data_script + \
              ' ' + input_soc_obj_fl + \
              ' ' + input_gene_GO_fl + \
              ' ' + log2fc_cutoff_val_final + \
              ' ' + input_prefix + \
              ' ' + open_GO_enrich_annot_final_dir
        print(cmd)
        subprocess.call(cmd, shell=True)








if __name__ == "__main__":
    main()


