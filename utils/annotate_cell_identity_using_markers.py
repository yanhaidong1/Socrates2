#!/usr/bin/env python

import argparse
import glob
import sys
import os
import subprocess

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

    parser.add_argument("-impute_gene_mtx", dest= 'impute_gene_acc_matrix', help = 'Provide the impute gene accessibility file.')

    parser.add_argument("-meta_fl", dest='meta_file', help='Provide a meta file aftering the clustering.')

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

    parser.add_argument("-prepared_gene_fl", dest = 'prepared_gene_file', help = 'Once we open only plot, we need to specify a prepared gene file used for the dotplot.')

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


    if args.open_markerUMAP is None:

        open_markerUMAP_final = 'no'

    else:
        if args.open_markerUMAP == 'yes':

            open_markerUMAP_final = 'yes'

            ##input meta file
            if args.meta_file is None:
                print('Cannot find meta file, please provide it')
                return
            else:
                try:
                    file = open(args.meta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the meta file!')
                    return

            ##check the gff file
            if args.impute_gene_acc_matrix is None:
                print('Cannot find imputed gene accessibility file, please provide it')
                return
            else:
                try:
                    file = open(args.impute_gene_acc_matrix, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the imputed gene accessibility file!')
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

                if args.gene_tn5_sparse_file is None:
                    print('Cannot find gene tn5 sparse file, please provide it')
                    return
                else:
                    try:
                        file = open(args.gene_tn5_sparse_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the gene tn5 sparse file!')
                        return

                if args.impute_gene_acc_matrix is None:
                    print('Cannot find imputed gene accessibility file, please provide it')
                    return
                else:
                    try:
                        file = open(args.impute_gene_acc_matrix, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the imputed gene accessibility file!')
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

                ##input meta file
                if args.meta_file is None:
                    print('Cannot find meta file, please provide it')
                    return
                else:
                    try:
                        file = open(args.meta_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the meta file!')
                        return

            else:

                if args.marker_file is None:
                    print('Cannot find marker gene file, please provide it')
                    return
                else:
                    try:
                        file = open(args.marker_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the marker gene file!')
                        return

                if args.prepared_gene_file is None:
                    print('Cannot find prepared gene file, please provide it')
                    return
                else:
                    try:
                        file = open(args.prepared_gene_file, 'r')  ##check if the file is not the right file
                    except IOError:
                        print('There was an error opening the prepared gene file!')
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

            if args.gene_accessibility_mtx is None:
                print('Cannot find gene accessibility file, please provide it')
                return
            else:
                try:
                    file = open(args.gene_accessibility_mtx, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the gene accessibility file!')
                    return


            ##input meta file
            if args.meta_file is None:
                print('Cannot find meta file, please provide it')
                return
            else:
                try:
                    file = open(args.meta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the meta file!')
                    return

            ##input svd file
            if args.svd_file is None:
                print('Cannot find svd file, please provide it')
                return
            else:
                try:
                    file = open(args.svd_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the svd file!')
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

        open_markerUMAP_final_dir = output_dir + '/open_markerUMAP_final_dir'
        if not os.path.exists(open_markerUMAP_final_dir):
            os.makedirs(open_markerUMAP_final_dir)


        vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_UMAP.R'
        ipt_meta_fl = args.meta_file
        ipt_impute_fl = args.impute_gene_acc_matrix
        ipt_marker_fl = args.marker_file

        cmd = 'Rscript ' + vis_R_script + \
              ' ' + ipt_meta_fl + \
              ' ' + ipt_impute_fl + \
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

            ipt_impute_fl = args.impute_gene_acc_matrix
            ipt_sparse_fl = args.gene_tn5_sparse_file
            ipt_meta_fl = args.meta_file

            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + ipt_impute_fl + \
                  ' ' + ipt_sparse_fl + \
                  ' ' + ipt_meta_fl + \
                  ' ' + target_cluster_nm_final + \
                  ' ' + open_vis_marker_dotplot_final_dir
            print(cmd)
            subprocess.call(cmd,shell=True)

            vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_dotplot_plotting.R'

            ipt_prepared_gene_fl = output_dir + '/open_vis_marker_dotplot_final_dir/opt_gene_zscore_accprop.rds'
            ipt_marker_fl = args.marker_file

            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + ipt_prepared_gene_fl + \
                  ' ' + ipt_marker_fl + \
                  ' ' + open_vis_marker_dotplot_final_dir + \
                  ' ' + dotplotwidth_val_final + \
                  ' ' + dotplotheight_val_final
            print(cmd)
            subprocess.call(cmd, shell=True)

        else:

            vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_markers_dotplot_plotting.R'

            ipt_prepared_gene_fl = output_dir + '/open_vis_marker_dotplot_final_dir/opt_gene_zscore_accprop.rds'
            ipt_marker_fl = args.marker_file

            cmd = 'Rscript ' + vis_R_script + \
                  ' ' + ipt_prepared_gene_fl + \
                  ' ' + ipt_marker_fl + \
                  ' ' + open_vis_marker_dotplot_final_dir + \
                  ' ' + dotplotwidth_val_final + \
                  ' ' + dotplotheight_val_final
            print(cmd)
            subprocess.call(cmd,shell=True)

    if open_aggregate_annotation_final == 'yes':

        print('Users choose to open the aggregating marker accessibility to visualize the markers')

        open_aggregate_annotation_final_dir = output_dir + '/open_aggregate_annotation_final_dir'
        if not os.path.exists(open_aggregate_annotation_final_dir):
            os.makedirs(open_aggregate_annotation_final_dir)

        vis_R_script = input_required_scripts_dir + '/utils_annotate_cell_identity_using_markers/visualizations_aggregate_markers_UMAP.R'

        ipt_gene_acc_mtx_rds_fl = args.gene_accessibility_mtx
        ipt_meta_fl = args.meta_file
        ipt_svd_fl = args.svd_file
        prefix = 'opt_annot'
        openAnnot = 'yes'
        ipt_marker_fl = args.marker_file

        cmd = 'Rscript ' + vis_R_script + \
              ' ' + ipt_gene_acc_mtx_rds_fl + \
              ' ' + ipt_meta_fl + \
              ' ' + ipt_marker_fl + \
              ' ' + ipt_svd_fl + \
              ' ' + prefix + \
              ' ' + openAnnot + \
              ' ' + open_aggregate_annotation_final_dir
        print(cmd)
        subprocess.call(cmd,shell=True)







if __name__ == "__main__":
    main()


