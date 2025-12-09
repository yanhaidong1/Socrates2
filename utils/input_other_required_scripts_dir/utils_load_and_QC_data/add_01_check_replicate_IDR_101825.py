#!/usr/bin/env python

import argparse
import glob
import sys
import os
import re
import subprocess
import pandas as pd


##this script is to check the replication for the IDR

input_rep1_peak_fl = sys.argv[1]
##output.findcell_obj.rds

input_rep2_peak_fl = sys.argv[2]
##output.findcell_obj.rds

input_prefix = sys.argv[3]

input_output_dir = sys.argv[4]

def perform_IDR (input_rep1_peak_fl,input_rep2_peak_fl,input_prefix,input_output_dir):

    cmd = 'idr ' + \
          ' --samples ' + input_rep1_peak_fl + ' ' + input_rep2_peak_fl + \
          ' --idr-threshold 0.05 ' + \
          ' --output-file ' + input_output_dir + '/opt_' + input_prefix + '_idr.txt' + \
          ' --plot ' + \
          ' --input-file-type ' + 'narrowPeak' + \
          ' --rank p.value 2> ' + input_output_dir + '/idr_summary.log'
          #' >' + input_output_dir + '/opt_' + input_prefix + '_log.txt' + ' 2>&1'
    print(cmd)
    subprocess.call(cmd,shell=True)

    output_file = input_output_dir + '/opt_' + input_prefix + '_idr.txt'

    if os.path.isfile(output_file):
        # 读取 narrowPeak 文件
        df = pd.read_csv(output_file, sep="\t", header=None)

        # 添加第4列：chr:start_end
        df[3] = df[0].astype(str) + ":" + df[1].astype(str) + "_" + df[2].astype(str)

        # 按染色体和坐标排序
        df = df.sort_values(by=[0, 1, 2])

        # 保存结果（不含行列名）
        df.to_csv(input_output_dir + '/opt_' + input_prefix + '_idr_narrowpeak.txt', sep="\t", index=False, header=False)

        #print(f"完成格式化: {sample_name}_peaks.narrowPeak")

    else:
        print("no file identified")


perform_IDR (input_rep1_peak_fl,input_rep2_peak_fl,input_prefix,input_output_dir)





