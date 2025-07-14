#!/usr/bin/env python

##this script we will make the peak sparse file for each organ with using the parallele

##updating 071425 debug do not replace chr if there is only chr

##this script is to call peaks for each lib
import re
import glob
import subprocess
import os



def make_sparse (input_peak_fl,input_tn5_bed_fl,input_fastSparsetn5_pl,
                 input_output_dir):

    prefix = 'opt'

    opt_sorted_peak_bed_dir = input_output_dir + '/opt_sorted_peak_bed_dir'
    if not os.path.exists(opt_sorted_peak_bed_dir):
        os.makedirs(opt_sorted_peak_bed_dir)

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'
    if not os.path.exists(opt_peak_sparse_dir):
        os.makedirs(opt_peak_sparse_dir)

    peak_bed_fl_path = input_peak_fl


    #peak_bed_fl_path = eachdir + '/' + dir_nm + '.unique500bpPeaks.bed'
    sorted_tn5_bed_fl_path = input_tn5_bed_fl

    cmd = 'sort -k1,1V -k2,2n ' + peak_bed_fl_path + ' > ' + opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)
    sorted_peak_bed_fl_path = opt_sorted_peak_bed_dir + '/' + prefix + '.unique500bpPeaks_sorted.bed'

    cmd = 'bedtools intersect -a ' +  sorted_tn5_bed_fl_path + \
          ' -b ' + sorted_peak_bed_fl_path + ' -wa -wb' + \
          ' -sorted | python ' + input_fastSparsetn5_pl + \
          ' - > ' + opt_peak_sparse_dir + '/opt_peak.sparse'
    print(cmd)
    subprocess.call(cmd,shell=True)


def sort_sparse (input_output_dir):

    opt_peak_sparse_dir = input_output_dir + '/opt_peak_sparse_dir'

    opt_peak_sparse_sorted_dir = input_output_dir + '/opt_peak_sparse_sorted_dir'
    if not os.path.exists(opt_peak_sparse_sorted_dir):
        os.makedirs(opt_peak_sparse_sorted_dir)

    allsparse_fl_list = glob.glob(opt_peak_sparse_dir + '/*')
    for eachfl in allsparse_fl_list:
        mt = re.match('.+/(.+)\.sparse',eachfl)
        flnm = mt.group(1)

        cmd = 'sort -k1,1V -k2,2n ' + eachfl + ' > ' + opt_peak_sparse_sorted_dir + '/' + flnm + '_sorted.sparse'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##remove the peak sparse file
        cmd = 'rm ' + eachfl
        print(cmd)
        subprocess.call(cmd,shell=True)


def prepare_peak_acc (input_peak_tn5_fl,input_final_peak_fl,input_output_dir):

    step01_prepare_peak_acc_dir = input_output_dir + '/opt_prepare_peak_acc_dir'
    if not os.path.exists(step01_prepare_peak_acc_dir):
        os.makedirs(step01_prepare_peak_acc_dir)



    store_final_peak_dic = {}
    with open (input_final_peak_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            final_line = col[0] + '_' + col[1] + '_' + col[2]
            store_final_peak_dic[final_line] = 1

    store_all_acr_dic = {}
    with open (input_peak_tn5_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            ##updating 071425
            peak_list = col[0].split('_')
            chrnm = peak_list[0]
            if re.match('chr\d+',chrnm):
                new_peak = col[0]
            else:
                new_peak = col[0].replace('chr','')
            acrloc = '\t'.join(new_peak.split('_'))
            store_all_acr_dic[acrloc] = new_peak

    store_acr_ID_dic = {}
    store_final_line_list = []
    count = 0
    for eachacrloc in store_all_acr_dic:

        acrloc_nospace = store_all_acr_dic[eachacrloc]

        count += 1
        acrID = 'scACR_' + str(count)

        if acrloc_nospace in store_final_peak_dic:

            final_line = eachacrloc + '\t' + acrID + '\t' + '1'
            store_final_line_list.append(final_line)

            store_acr_ID_dic[acrloc_nospace] = acrID

    with open (step01_prepare_peak_acc_dir + '/opt_peak.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    tn5_fl = open (input_peak_tn5_fl,'r')
    opt_fl = open (step01_prepare_peak_acc_dir + '/opt_accessibility.txt','w')

    for eachline in tn5_fl:

        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        acrlocnospace = col[0]
        acrlocnospace = acrlocnospace.replace('chr','')

        if acrlocnospace in store_final_peak_dic:

            acrID = store_acr_ID_dic[acrlocnospace]

            final_line = acrID + '\t' + col[1] + '\t' + col[2]
            opt_fl.write(final_line + '\n')

    tn5_fl.close()
    opt_fl.close()




