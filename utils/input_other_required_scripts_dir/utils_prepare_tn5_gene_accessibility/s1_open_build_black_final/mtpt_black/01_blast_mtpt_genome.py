#!/usr/bin/env python

##this script we will blast mtpt to the genome
import re
import glob
import sys
import subprocess
import os
from Bio import SeqIO

input_genome_fl = sys.argv[1]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/Osativa323v7.fa

input_output_dir = sys.argv[2]

input_organelle_str = sys.argv[3]

def make_blast (input_genome_fl,input_output_dir,input_organelle_str):

    store_mtpt_genome_seq_dic = {}
    for seq_record in SeqIO.parse(input_genome_fl,'fasta'):
        if seq_record.id in input_organelle_str.split(','):
            store_mtpt_genome_seq_dic[seq_record.id] = str(seq_record.seq)

    with open (input_output_dir + '/chrptmt.fa','w+') as opt:
        for eachchr in store_mtpt_genome_seq_dic:
            opt.write('>' + eachchr + '\n' + store_mtpt_genome_seq_dic[eachchr] + '\n')

    ##make a blast
    cmd = 'makeblastdb -in ' + input_genome_fl + ' –dbtype nucl –out ' + input_output_dir + '/db'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'blastn ' +  '-db ' + input_output_dir + '/db' + ' -query ' + input_output_dir + '/chrptmt.fa' + \
          ' -out ' + input_output_dir + '/opt_blast.txt' + \
          ' -evalue 0.00001 -outfmt 6'
    subprocess.call(cmd,shell=True)

make_blast (input_genome_fl,input_output_dir,input_organelle_str)
