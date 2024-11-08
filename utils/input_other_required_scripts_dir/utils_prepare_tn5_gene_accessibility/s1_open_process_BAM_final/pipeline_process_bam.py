#!/usr/bin/env python

##this script is to process bam files from the cell ranger output

import subprocess


def process_bam (input_other_required_scripts_dir,input_bam_fl,input_output_dir,
                 input_core_num,input_qual_val,remove_temp_file):


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

    cmd = 'java -jar $EBROOTPICARD/picard.jar' + \
          ' MarkDuplicates' + \
          ' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000' + \
          ' REMOVE_DUPLICATES=true' + \
          ' I=' + input_output_dir + '/temp_mapped_sorted.bam' + \
          ' O=' + input_output_dir + '/temp_mapped_rmpcr.bam'+ \
          ' METRICS_FILE=' + input_output_dir + '/dups.txt' + \
          ' BARCODE_TAG=CB' + \
          ' ASSUME_SORT_ORDER=coordinate' + \
          ' USE_JDK_DEFLATER=true' + \
          ' USE_JDK_INFLATER=true'
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

    ##filter bam MQ > qual and properly paired
    cmd = 'samtools view -@ ' + input_core_num + ' -f 3 -bhq ' + input_qual_val + ' ' + input_output_dir + '/temp_mapped_rmpcr.bam > ' + input_output_dir + '/temp_mq' + input_qual_val + '_rmpcr.bam'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##filter barcodes
    cmd = 'perl ' + input_target_required_scripts_dir + '/fixBC.pl ' + input_output_dir + '/temp_mq' + input_qual_val + '_rmpcr.bam | samtools view -bhS - > ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##make Tn5 bed files
    cmd = 'perl ' + input_target_required_scripts_dir + '/makeTn5bed.pl ' + input_output_dir + '/temp_fixBC_mq' + input_qual_val + '_rmpcr.bam | sort -k1,1 -k2,2n - > ' + input_output_dir + '/opt_tn5_mq' + input_qual_val + '.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##remove the temp files
    if remove_temp_file == 'yes':

        cmd = 'rm ' + input_output_dir + '/*bam'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'rm ' + input_output_dir + '/dups.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)



