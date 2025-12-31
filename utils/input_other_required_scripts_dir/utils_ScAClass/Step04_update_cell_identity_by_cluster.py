#!/usr/bin/env python

##updating 101325 we will update the cluster for all the examined method
##this script is to update the cell identity by cluster information

import re
import glob
import sys
import subprocess
import os


input_clustering_fl = sys.argv[1]
input_Rscript_fl = sys.argv[2]
input_cell_prediction_fl = sys.argv[3]
input_output_dir = sys.argv[4]

def update_prediction_by_cluster (input_clustering_fl,input_Rscript_fl,input_cell_prediction_fl,input_output_dir):


    cmd = 'Rscript ' + input_Rscript_fl + \
          ' ' + input_clustering_fl + \
          ' ' + input_output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

    ipt_cluster_meta_fl = input_output_dir + '/temp_meta_fl.txt'


    ##eachPC + '\t' + eachmethod + '\t' + eachcell + '\t' + final_annot
    store_cluster_cell_list_dic = {}
    count = 0
    with open (ipt_cluster_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cellnm = col[0]
            cluster = col[-1]
            count += 1
            if count != 1:
                if cluster in store_cluster_cell_list_dic:
                    store_cluster_cell_list_dic[cluster].append(cellnm)
                else:
                    store_cluster_cell_list_dic[cluster] = []
                    store_cluster_cell_list_dic[cluster].append(cellnm)

    ##we will first check the methods output in this file
    all_method_dic = {}
    count = 0
    with open(input_cell_prediction_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                method = col[1]
                all_method_dic[method] = 1

    for eachmethod in all_method_dic:


        store_cell_annot_dic = {}
        with open (input_cell_prediction_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                cellnm = col[2]
                annot = col[3]

                method = col[1]
                if method == eachmethod:

                    ##updating 122125
                    ##we will not consider the Unknown
                    if annot != 'Unknown':

                        store_cell_annot_dic[cellnm] = annot

        store_final_line_list = []
        store_final_line_eachcell_list = []
        for eachcluster in store_cluster_cell_list_dic:
            cluster_cellnm_list = store_cluster_cell_list_dic[eachcluster]
            total_cell_num = len(cluster_cellnm_list)

            store_celltype_prop_dic = {}
            ##we will check the proportion of each cell type
            store_celltype_numcell_dic = {}
            for eachcell in store_cell_annot_dic:

                if eachcell in cluster_cellnm_list:
                    celltype = store_cell_annot_dic[eachcell]

                    if celltype in store_celltype_numcell_dic:
                        store_celltype_numcell_dic[celltype] += 1
                    else:
                        store_celltype_numcell_dic[celltype] = 1

            for eachcelltype in store_celltype_numcell_dic:
                cellnum = store_celltype_numcell_dic[eachcelltype]
                prop = cellnum/total_cell_num
                store_celltype_prop_dic[eachcelltype] = prop

            ##updating 122125
            if len(list(store_celltype_prop_dic.keys())) != 0:

                ##find the max prop for the eachcelltype
                max_celltype = max(store_celltype_prop_dic, key=store_celltype_prop_dic.get)

                max_celltype_prop = store_celltype_prop_dic[max_celltype]

                final_line = eachcluster + '\t' + max_celltype + '\t' + str(max_celltype_prop)
                store_final_line_list.append(final_line)

                for eachcell in cluster_cellnm_list:
                    final_line = eachcluster + '\t' + eachcell + '\t' + max_celltype + '\t' + str(max_celltype_prop)
                    store_final_line_eachcell_list.append(final_line)

        with open (input_output_dir + '/opt_' + eachmethod + '_final_celltype_annotation_per_cluster.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        with open(input_output_dir + '/opt_final_' + eachmethod + '_cell_annotation_per_cluster.txt', 'w+') as opt:
            for eachline in store_final_line_eachcell_list:
                opt.write(eachline + '\n')

        ##updating 121324 we will update the cell identity in the original cluster file
        store_cell_annot_not_by_cluster_dic = {}
        count = 0
        with open (input_output_dir + '/opt_final_cell_annotation.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    cellID = col[2]
                    annot = col[3]
                    method = col[1]

                    if method == eachmethod:

                        store_cell_annot_not_by_cluster_dic[cellID] = annot

        store_cell_annot_by_cluster_dic = {}
        with open (input_output_dir + '/opt_final_' + eachmethod + '_cell_annotation_per_cluster.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                cellID = col[1]
                annot = col[2]
                store_cell_annot_by_cluster_dic[cellID] = annot

        store_final_line_list = []
        count = 0
        with open (ipt_cluster_meta_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    cellID = col[0]
                    if cellID in store_cell_annot_not_by_cluster_dic:
                        ScAClass = store_cell_annot_not_by_cluster_dic[cellID]
                    else:
                        ScAClass = 'Unknown'

                    if cellID in store_cell_annot_by_cluster_dic:
                        ScAClass_clust = store_cell_annot_by_cluster_dic[cellID]
                    else:
                        ScAClass_clust = 'Unknown'

                    final_line = eachline + '\t' + ScAClass + '\t' + ScAClass_clust
                    store_final_line_list.append(final_line)
                else:
                    final_line = eachline + '\t' + 'ScATACtor' + '\t' + 'ScATACtor_clust'
                    store_final_line_list.append(final_line)

        with open (input_output_dir + '/opt_' + eachmethod + '_meta_ScATACtor_predicted_cell_identity.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

update_prediction_by_cluster (input_clustering_fl,input_Rscript_fl,input_cell_prediction_fl,input_output_dir)










