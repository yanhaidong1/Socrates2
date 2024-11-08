#!/usr/bin/env python

##this script is to update the cell identity by cluster information

import re
import glob
import sys
import subprocess
import os

input_clustering_fl = sys.argv[1]
input_cell_prediction_fl = sys.argv[2]
input_output_dir = sys.argv[3]

def update_prediction_by_cluster (input_clustering_fl,input_cell_prediction_fl,input_output_dir):

    ##eachPC + '\t' + eachmethod + '\t' + eachcell + '\t' + final_annot
    store_cluster_cell_list_dic = {}
    with open (input_clustering_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cellnm = col[0]
            cluster = col[1]
            if cluster in store_cluster_cell_list_dic:
                store_cluster_cell_list_dic[cluster].append(cellnm)
            else:
                store_cluster_cell_list_dic[cluster] = []
                store_cluster_cell_list_dic[cluster].append(cellnm)

    store_cell_annot_dic = {}
    with open (input_cell_prediction_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cellnm = col[2]
            annot = col[3]
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

        ##find the max prop for the eachcelltype
        max_celltype = max(store_celltype_prop_dic, key=store_celltype_prop_dic.get)

        max_celltype_prop = store_celltype_prop_dic[max_celltype]

        final_line = eachcluster + '\t' + max_celltype + '\t' + str(max_celltype_prop)
        store_final_line_list.append(final_line)

        for eachcell in cluster_cellnm_list:
            final_line = eachcluster + '\t' + eachcell + '\t' + max_celltype + '\t' + str(max_celltype_prop)
            store_final_line_eachcell_list.append(final_line)

    with open (input_output_dir + '/opt_final_celltype_annotation_per_cluster.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    with open(input_output_dir + '/opt_final_cell_annotation_per_cluster.txt', 'w+') as opt:
        for eachline in store_final_line_eachcell_list:
            opt.write(eachline + '\n')


update_prediction_by_cluster (input_clustering_fl,input_cell_prediction_fl,input_output_dir)










