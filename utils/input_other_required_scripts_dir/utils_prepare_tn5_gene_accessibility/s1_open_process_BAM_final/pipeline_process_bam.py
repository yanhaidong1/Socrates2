#!/usr/bin/env python

##updating 010526 we will set an option to remove PCR dup in another way
##updating 010226 we will check if the bam includes the RG
##updating 021825 we will modify the bam to remove the -1
##updating 012325 add a new function to remove the black
##this script is to process bam files from the cell ranger output

import subprocess
import re
import os
import pysam



def process_bam (input_other_required_scripts_dir,input_bam_fl,input_output_dir,
                 input_core_num,input_qual_val,remove_temp_file,use_picard_tool_final):


    input_target_required_scripts_dir = input_other_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_process_BAM_final'


    ##retain only mapped reads
    cmd = 'samtools view -@ ' + str(input_core_num) + ' -bhq 1 ' + input_bam_fl + ' > ' + input_output_dir + '/temp_mapped.bam'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##sort the bam
    cmd = 'samtools' + \
          ' sort ' + input_output_dir + '/temp_mapped.bam'+ \
          ' -@ ' + input_core_num + \
          ' -o ' + input_output_dir + '/temp_mapped_sorted.bam' + \
          ' > ' + input_output_dir + '/temp_mapped_sorted.bam'
    print(cmd)
    subprocess.call(cmd, shell=True)

    input_bam = input_output_dir + '/temp_mapped_sorted.bam'
    output_bam = input_output_dir + '/temp_mapped_sorted_modi.bam'

    with pysam.AlignmentFile(input_bam, "rb") as infile, \
            pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:

        for read in infile:
            if read.has_tag("CB"):
                barcode = read.get_tag("CB").rstrip("-1")
                read.set_tag("CB", barcode, value_type='Z', replace=True)
            outfile.write(read)

    ##cmd = 'samtools view -h ' + input_output_dir + '/temp_mapped_sorted.bam' + ' | awk \'{if ($0 ~ /^@/) {print; next} for (i=1; i<=NF; i++) {if ($i ~ /^CB:Z:/) sub(/-1$/, "", $i)} print}\' ' \
    ##      + ' | samtools view -bh -o ' +input_output_dir + '/temp_mapped_sorted_modi.bam'
    ##print(cmd)
    ##ubprocess.call(cmd,shell=True)

    ##updating 010226
    cmd = 'samtools view -H ' + output_bam + ' | grep \'^@RG\' > ' + input_output_dir + '/temp_header.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    if os.path.getsize( input_output_dir + '/temp_header.txt') == 0:
        print("The BAM file does not include RG information, and this pipeline will add the RG to the BAM")

        cmd = 'samtools addreplacerg -r ID:rg1 -r SM:sample1 -r PL:ILLUMINA -o ' + input_output_dir + '/temp_temp_mapped_sorted_modi_addRG.bam' + ' ' + input_output_dir + '/temp_mapped_sorted_modi.bam'
        print(cmd)
        subprocess.call(cmd,shell=True)

        #final_bam = input_output_dir + '/temp_temp_mapped_sorted_modi_addRG.bam'

        cmd = 'samtools' + \
              ' sort ' + input_output_dir + '/temp_temp_mapped_sorted_modi_addRG.bam' + \
              ' -@ ' + input_core_num + \
              ' -o ' + input_output_dir + '/temp_temp_mapped_sorted_modi_addRG_sorted.bam' + \
              ' > ' + input_output_dir + '/temp_temp_mapped_sorted_modi_addRG_sorted.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)

        final_bam = input_output_dir + '/temp_temp_mapped_sorted_modi_addRG_sorted.bam'



    else:
        final_bam = input_output_dir + '/temp_mapped_sorted_modi.bam'


    if use_picard_tool_final == 'yes':

        cmd = 'picard' + \
              ' MarkDuplicates' + \
              ' --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000' + \
              ' --REMOVE_DUPLICATES true' + \
              ' -I ' + final_bam + \
              ' -O ' + input_output_dir + '/temp_mapped_rmpcr.bam'+ \
              ' --METRICS_FILE ' + input_output_dir + '/dups.txt' + \
              ' --BARCODE_TAG CB' + \
              ' --ASSUME_SORT_ORDER coordinate' + \
              ' --USE_JDK_DEFLATER true' + \
              ' --USE_JDK_INFLATER true'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##run picard to remove the PCR duplicated reads
        #cmd = 'java -jar $EBROOTPICARD/picard.jar' + \
        #      ' MarkDuplicates' + \
        #      ' -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000' + \
        #      ' -REMOVE_DUPLICATES true' + \
        #      ' -I ' + input_output_dir + '/temp_mapped.bam' + \
        #      ' -O ' + input_output_dir + '/temp_mapped_rmpcr.bam'+ \
        #      ' -M ' + input_output_dir + '/dups.txt' + \
        #      ' -BARCODE_TAG CB' + \
        #      ' -ASSUME_SORT_ORDER coordinate'
              #' VALIDATION_STRINGENCY=LENIENT'
              #' ASSUME_SORTED=true' + \
              ##MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 is not valide
        #subprocess.call(cmd, shell=True)

    else:

        print('Users do not want to use picard to remove PCR duplicates')
        ipt_python_script =  input_other_required_scripts_dir + '/utils_prepare_tn5_gene_accessibility/s1_open_process_BAM_final/scATACseq_cb_pcr_dedup.py'

        cmd = 'python ' + ipt_python_script + \
              ' -o ' + input_output_dir + '/temp_mapped_rmpcr.bam' + \
              ' -i ' + final_bam + \
              ' -p ' + input_core_num
        print(cmd)
        subprocess.call(cmd,shell=True)






    ##filter bam MQ > qual and properly paired
    cmd = 'samtools view -@ ' + input_core_num + ' -f 3 -bhq ' + input_qual_val + ' ' + input_output_dir + '/temp_mapped_rmpcr.bam > ' + input_output_dir + '/temp_mq' + input_qual_val + '_rmpcr.bam'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##filter barcodes
    cmd = 'perl ' + input_target_required_scripts_dir + '/fixBC.pl ' + input_output_dir + '/temp_mq' + input_qual_val + '_rmpcr.bam | samtools view -bhS - > ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##updating 042325 use the python script other than the pl
    #cmd = 'python ' + input_target_required_scripts_dir + '/fixBC.py ' + input_output_dir + '/temp_mq' + input_qual_val + '_rmpcr.bam | samtools view -bhS - > ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam'
    #print(cmd)
    #subprocess.call(cmd,shell=True)


    ##make Tn5 bed files
    current_directory = os.getcwd()
    cmd = 'cp ' + input_target_required_scripts_dir + '/perlSAM.pm ' + current_directory
    print(cmd)
    subprocess.call(cmd,shell=True)


    cmd = 'perl ' + input_target_required_scripts_dir + '/makeTn5bed.pl ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam | sort -k1,1 -k2,2n - > ' + input_output_dir + '/opt_tn5_mq' + input_qual_val + '.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##udpating 042325 use the python script other than the pl
    #cmd = 'python ' + input_target_required_scripts_dir + '/makeTn5bed_transfer.py ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam | sort -k1,1 -k2,2n - > ' + input_output_dir + '/opt_tn5_mq' + input_qual_val + '.bed'
    #print(cmd)
    #subprocess.call(cmd, shell=True)


    ##remove the temp files
    if remove_temp_file == 'yes':

        cmd = 'rm ' + input_output_dir + '/*bam'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'rm ' + input_output_dir + '/dups.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        if os.path.exists(current_directory + '/perlSAM.pm'):
            cmd = 'rm ' + current_directory + '/perlSAM.pm'
            print(cmd)
            subprocess.call(cmd, shell=True)



def remove_black (ipt_tn5_fl,input_black_list_fl,opt_dir,remove_temp_file):

    if re.match('.+/.+',ipt_tn5_fl):
        mt = re.match('.+/(.+)\.bed',ipt_tn5_fl)
        flnm = mt.group(1)
    else:
        mt = re.match('(.+)\.bed',ipt_tn5_fl)
        flnm = mt.group(1)

    ##remove the black list
    cmd = 'sort -k1,1V -k2,2n ' + input_black_list_fl + ' > ' + opt_dir + '/temp_black_list_sorted.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##uniq
    cmd = 'uniq ' + ipt_tn5_fl + ' > ' + opt_dir + '/temp_uniq_tn5.bed'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##sort
    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_uniq_tn5.bed' + ' > ' + opt_dir + '/temp_uniq_tn5_sorted.bed'
    print(cmd)
    subprocess.call(cmd, shell=True)


    ##intersect with the bed
    cmd = 'bedtools intersect -a ' + opt_dir + '/temp_uniq_tn5_sorted.bed' + ' -b ' + opt_dir + '/temp_black_list_sorted.bed' + \
          ' -wa -wb ' + ' > ' + opt_dir + '/temp_intersect.txt'
    subprocess.call(cmd, shell=True)

    store_black_region_dic = {}
    with open(opt_dir + '/temp_intersect.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            black_loc = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3] + '_' + col[4]
            store_black_region_dic[black_loc] = 1

    bed_file = open(opt_dir + '/temp_uniq_tn5_sorted.bed', 'r')
    bed_file_noblack_fl = open(opt_dir + '/' + flnm + '_rmblack.bed', 'w')

    for eachline in bed_file:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()
        loc = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3] + '_' + col[4]
        if loc not in store_black_region_dic:
            bed_file_noblack_fl.write(eachline + '\n')

    bed_file.close()
    bed_file_noblack_fl.close()

    ##remove the temp files
    if remove_temp_file == 'yes':

        cmd = 'rm ' + opt_dir + '/temp*'
        print(cmd)
        subprocess.call(cmd,shell=True)







