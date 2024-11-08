#!/usr/bin/env python

##updating 061923 we will set the user provided meta
##updating 112821 do the seperate annotation
##updating 110221 we need to plot the UMAP from the prediction
##this script is to do the clustering annotation of each cluster
import re
import glob
import sys
import subprocess
import os



input_store_clusters_dir = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/06_1_seperate_clustering_112421

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/04_1_subclustering_071821/output_dir_NMF_Socrates_102421/store_split_cluster_dir

input_store_norm_clusters_dir = sys.argv[2]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/06_2_seperate_norm_acc_visualize_markers_112421

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/04_2_norm_acc_visualize_markers_subclust_102221/output_dir_Socrate_nvar_NMF_102221/

input_marker_dir = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/03_3_norm_acc_visualize_markers_071421/input_marker_dir

input_annotate_script = sys.argv[4]
##annotate_cells.v3_noflt_markers.R

input_plot_UMAP_script = sys.argv[5]
##01_plot_UMAP_clustering_annotation_110221.R

input_config_fl = sys.argv[6]

input_output_dir = sys.argv[7] ##output_dir



store_target_parameter_dic = {}
with open (input_config_fl,'r') as ipt:
    for eachline in ipt:
        eachline = eachline.strip('\n')
        if not eachline.startswith('#'):
            col = eachline.strip().split('=')
            store_target_parameter_dic[col[0]] = col[1]


#clustannot = store_target_parameter_dic['clustannot']
#target_cluster_string = store_target_parameter_dic['target_cluster_string']
#target_cluster_list = target_cluster_string.split(',')
openAnnot = store_target_parameter_dic['openAnnot'] ##yes or no
##updating 110221 for the annotation plotting
target_cluster_annot_str = store_target_parameter_dic['target_cluster_annot_str'] ##it is cluster_annotation_smooth,cluster_annotation_knn,cluster_annotation_glmnet
target_cluster_str_annot_list = target_cluster_annot_str.split(',')

target_tissue_str = store_target_parameter_dic['target_tissue_str']
target_tissue_str_list = target_tissue_str.split(',')

##updating 061923
target_parameter_setting_str_forcluster0601 = store_target_parameter_dic['target_parameter_setting_str_forcluster0601'] ##win10k_PCs30_NMF_dir
target_parameter_setting_str_forvis0602 = store_target_parameter_dic['target_parameter_setting_str_forvis0602']

open_user_provided_meta_fl = store_target_parameter_dic['open_user_provided_meta_fl']



def pipeline_analyze (input_store_clusters_dir,input_store_norm_clusters_dir,input_marker_dir,input_output_dir,
                      input_annotate_script,openAnnot,target_tissue_str_list,target_parameter_setting_str_forcluster0601,target_parameter_setting_str_forvis0602,open_user_provided_meta_fl):



    cluster_dir_list = glob.glob(input_store_clusters_dir + '/*')

    for eachcluster_dir in cluster_dir_list:
        mt = re.match('.+/(.+)',eachcluster_dir)
        tissue_nm = mt.group(1)

        print('tissue name is ' + tissue_nm)
        print(target_tissue_str_list)

        if tissue_nm in target_tissue_str_list:

            print('target tissue name is ' + tissue_nm)

            norm_data_fl = input_store_norm_clusters_dir + '/' + tissue_nm + '/output_dir_' + target_parameter_setting_str_forvis0602 + '/step01_firstround_create_obj_dir/all.normalizedActivity.sparse'
            #norm_data_fl = input_store_norm_clusters_dir + '/' + tissue_nm + '/step01_firstround_create_obj_dir/all.normalizedActivity.sparse'

            input_marker_fl = ''
            marker_fl_list = glob.glob(input_marker_dir + '/*')
            print(marker_fl_list)

            for eachmarker_fl in marker_fl_list:
                mt = re.match('.+/(.+)', eachmarker_fl)
                flnm = mt.group(1)
                print(flnm)
                mt = re.match('markers_(.+)_.+\.txt', flnm)
                marker_tissuenm = mt.group(1)

                #if 'Lseed' in tissue_nm or 'Eseed' in tissue_nm:
                #    new_tissue_nm = 'seed'
                #else:
                new_tissue_nm = tissue_nm

                if marker_tissuenm == new_tissue_nm:
                    input_marker_fl = eachmarker_fl

                    print('target marker file is')
                    print(input_marker_fl)

            store_tissue_dir = input_output_dir + '/' + tissue_nm
            if not os.path.exists(store_tissue_dir):
                os.makedirs(store_tissue_dir)

            ##find the meta file
            fl_list = glob.glob(eachcluster_dir + '/output_dir_' + target_parameter_setting_str_forcluster0601 + '/*')
            input_meta_fl = ''
            input_svd_fl = ''
            for eachfl in fl_list:
                mt = re.match('.+/(.+)', eachfl)
                flnm = mt.group(1)

                if re.match('.+res.+_metadata\.txt',flnm):
                    input_meta_fl = eachfl

                if 'reduced_dimensions.txt' in flnm:
                    input_svd_fl = eachfl

            if open_user_provided_meta_fl == 'no':
                input_final_meta_fl = input_meta_fl
            else:
                input_final_meta_fl = open_user_provided_meta_fl

            ##run the annotation script
            cmd = 'Rscript ' + input_annotate_script + \
                  ' ' + norm_data_fl + \
                  ' ' + input_final_meta_fl + \
                  ' ' + input_marker_fl + \
                  ' ' + input_svd_fl + \
                  ' ' + tissue_nm + \
                  ' ' + openAnnot + \
                  ' ' + store_tissue_dir
            print(cmd)
            subprocess.call(cmd,shell=True)


##updating 1110221
def pipeline_analysis_annot_UMAP (input_plot_UMAP_script,input_output_dir,target_cluster_str_annot_list):

    cluster_dir_list = glob.glob(input_output_dir + '/*')
    for eachclusterdir in cluster_dir_list:
        mt = re.match('.+/(.+)',eachclusterdir)
        tissue_nm = mt.group(1)

        ipt_all_meta_fl = eachclusterdir + '/' + tissue_nm + '.meta.txt'

        for eachtarget_col in target_cluster_str_annot_list:

            cmd = 'Rscript ' + input_plot_UMAP_script + \
                  ' ' + ipt_all_meta_fl + \
                  ' ' + eachtarget_col + \
                  ' ' + eachtarget_col + \
                  ' ' + eachtarget_col + \
                  ' ' + eachclusterdir
            print(cmd)
            subprocess.call(cmd,shell=True)


##open it when necessary
pipeline_analyze (input_store_clusters_dir,input_store_norm_clusters_dir,input_marker_dir,input_output_dir,
                      input_annotate_script,openAnnot,target_tissue_str_list,target_parameter_setting_str_forcluster0601,target_parameter_setting_str_forvis0602,open_user_provided_meta_fl)

pipeline_analysis_annot_UMAP (input_plot_UMAP_script,input_output_dir,target_cluster_str_annot_list)





