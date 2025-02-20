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

    parser.add_argument("-soc_obj", dest='soc_object_fl',
                        help='Provide an object obtained from the Calculate gene accessibility step.')


    ##Optional parameters
    parser.add_argument("-prop_threshold", dest='proportion_threshold',
                        help='selecting genes passing minimum proportion threshold.'
                             'Default: 0.01')

    parser.add_argument("-top_gene_num", dest = 'top_gene_number',
                        help='Decide the number of genes will be plotted per cluster.'
                             'Default: 5')

    parser.add_argument("-normT", dest = 'normalizeT',help = 'Determine the ranking method (mean.dif/prop.dif/adj.dif) for ordering genes based on top-to-bottom chromatin accessibility across different clusters.'
                                                             'Default: mean.dif.')

    parser.add_argument('-prefix', dest='prefix_name', help='Specify a prefix name for the final output file.'
                                                            'Default: output')


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
        print('Cannot find the soc object file, please provide it')
        return
    else:
        try:
            file = open(args.soc_object_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the soc object file!')
            return


    store_final_parameter_line_list = []

    if args.proportion_threshold is not None:
        proportion_threshold_final = args.proportion_threshold
    else:
        proportion_threshold_final = '0.01'

    store_final_parameter_line_list.append('threshold_prop <- ' + proportion_threshold_final)

    if args.top_gene_number is not None:
        top_gene_number_final = args.top_gene_number
    else:
        top_gene_number_final = '5'

    store_final_parameter_line_list.append('top_gene <- ' + top_gene_number_final)

    if args.normalizeT is not None:
        normalizeT_final = args.normalizeT
    else:
        normalizeT_final = 'mean.dif'

    store_final_parameter_line_list.append('normT_type <- ' + '\'' + normalizeT_final + '\'')

    if args.prefix_name is not None:
        prefix_name_final = args.prefix_name
    else:
        prefix_name_final = 'output'

    store_final_parameter_line_list.append('input_prefix <- ' + '\'' + prefix_name_final + '\'')


    ##we will firstly build up the parameter setting file
    with open(output_dir + '/temp_defined_parameters.config', 'w+') as opt:
        for eachline in store_final_parameter_line_list:
            opt.write(eachline + '\n')


    ipt_R_script = input_required_scripts_dir + '/utils_diff_markers_finding/rank_diff_genes.R'
    ipt_soc_obj_fl = args.soc_object_fl

    cmd = 'Rscript ' + ipt_R_script + \
          ' ' + output_dir + '/temp_defined_parameters.config' + \
          ' ' + ipt_soc_obj_fl + \
          ' ' + output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)








if __name__ == "__main__":
    main()


