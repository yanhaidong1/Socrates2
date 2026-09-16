#!/usr/bin/env python


import subprocess


##this script is to smooth gene accessibility


def smooth_gene_accessibility (input_required_scripts_dir,
                               input_soc_obj_fl,
                               input_output_dir,input_core_num,
                               target_cluster,input_prefix):

    gene_smooth_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s4_open_smooth_gene_accessibility_final/smooth_genes_acc_100724.R'
    function_script = input_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s4_open_smooth_gene_accessibility_final/subfunction_smooth_genes_acc_112723.R'

    cmd = 'Rscript ' + gene_smooth_script + \
          ' ' + input_soc_obj_fl + \
          ' ' + input_core_num + \
          ' ' + target_cluster + \
          ' ' + input_output_dir + \
          ' ' + function_script + \
          ' ' + input_prefix
    print(cmd)
    subprocess.call(cmd, shell=True)
