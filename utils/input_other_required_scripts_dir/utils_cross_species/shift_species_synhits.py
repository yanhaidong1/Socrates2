#!/usr/bin/env python


import argparse
import glob
import sys
import subprocess
import os
import re

input_synhit_fl = sys.argv[1]

input_output_dir = sys.argv[2]

ipt_prefix = sys.argv[3]

def shift_species (input_synhit_fl,input_output_dir,ipt_prefix):

    store_final_line_list = []
    count = 0
    with open (input_synhit_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:

                spe1_list = []
                for i in range(8,16):
                    spe1_list.append(col[i])
                spe1_str = '\t'.join(spe1_list)

                spe2_list = []
                for i in range(8):
                    spe2_list.append(col[i])

                spe2_str = '\t'.join(spe2_list)

                other_list = []
                for i in range(16,20):
                    other_list.append(col[i])

                other_str = '\t'.join(other_list)

                if len(col) == 21:

                    mt = re.match('(.+)_vs_(.+):(.+)', col[-1])
                    regID = mt.group(2) + '_vs_' + mt.group(1) + ':' + mt.group(3)
                    final_line = spe1_str + '\t' + spe2_str + '\t' + other_str + '\t' + regID
                    store_final_line_list.append(final_line)

                else:
                    ##for the col[-2]
                    mt = re.match('(.+)_vs_(.+):(.+)',col[-2])
                    regID = mt.group(2) + '_vs_' + mt.group(1) + ':' + mt.group(3)

                    mt = re.match('(.+)_vs_(.+):(.+)', col[-1])
                    blkID = mt.group(2) + '_vs_' + mt.group(1) + ':' + mt.group(3)

                    final_line = spe1_str + '\t' + spe2_str + '\t' + other_str + '\t' + regID + '\t' + blkID
                    store_final_line_list.append(final_line)

            else:
                store_final_line_list.append(eachline)

    with open (input_output_dir + '/' + ipt_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



shift_species (input_synhit_fl,input_output_dir,ipt_prefix)



