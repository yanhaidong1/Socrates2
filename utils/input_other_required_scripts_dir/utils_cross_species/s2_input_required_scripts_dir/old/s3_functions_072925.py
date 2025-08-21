#!/usr/bin/env python

##updating 053124 we will add the cell type information

import argparse
import os
import re



def modify_blasted_file (ipt_ori_blast_fl, ipt_spe1_ct_ACR_fl,ipt_spe2_syn_fl,opt_dir):

    ##add the ACR ID
    #store_spe1_ACR_loc_ACRID_dic = {}
    #ACR_ID = 0
    #with open (ipt_spe1_ACR_fl,'r') as ipt:
    #    for eachline in ipt:
    #        eachline = eachline.strip('\n')
    #        col = eachline.strip().split()
    #        ACRloc = col[0] + '_' + col[1] + '_' + col[2]
    #        ACR_ID += 1

    #        ACRID = 'scACR_' + str(ACR_ID)
    #        store_spe1_ACR_loc_ACRID_dic[ACRloc] = ACRID


    store_spe1_ACR_loc_ACRID_dic = {}
    store_spe1_ACR_loc_CT_dic = {}
    with open(ipt_spe1_ct_ACR_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            ACRloc = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('(.+);(.+)',col[3])
            ACRID = mt.group(1)
            celltype_str = mt.group(2)
            store_spe1_ACR_loc_ACRID_dic[ACRloc] = ACRID
            store_spe1_ACR_loc_CT_dic[ACRloc] = celltype_str


    store_spe2_syntenic_ID_loc_dic = {}
    with open (ipt_spe2_syn_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            syntenic_ID = col[-1]
            syntenic_loc = col[0] + ':' + col[1] + '-' + col[2]

            store_spe2_syntenic_ID_loc_dic[syntenic_ID] = syntenic_loc

    store_final_line_list = []
    with open (ipt_ori_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            ##updating 072925
            ##check if the blasted file contains '::'
            #if re.match('')

            syntenic_ID = col[1]
            spe2_syntenic_region = store_spe2_syntenic_ID_loc_dic[syntenic_ID]

            ACR_loc = col[0]
            ACR_loc_col = ACR_loc.split('_')
            final_ACR_loc = ACR_loc_col[0] + ':' + ACR_loc_col[1] + '-' + ACR_loc_col[2]
            ACRID = store_spe1_ACR_loc_ACRID_dic[ACR_loc]
            celltype = store_spe1_ACR_loc_CT_dic[ACR_loc]

            first_part = ACRID + ';' + celltype + ';' + syntenic_ID + '::' + final_ACR_loc + '\t' + \
                         syntenic_ID + '::' + spe2_syntenic_region

            second_part_list = []
            for i in range(2,len(col)):
                second_part_list.append(col[i])

            second_part_str = '\t'.join(second_part_list)

            final_line = first_part + '\t' + second_part_str
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_blasted_modif_format.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')






def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass


def process_blast_file(blast_file, evalue_threshold, output_file):
    passing_blast_lines = []
    passing_bed_lines = []
    passing_ref_bed_lines = []
    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split("\t")
            if int(columns[3]) > 20 and float(columns[10]) < evalue_threshold:
                # Extract the 'syntenic_regions####' part from the qseqid and sseqid
                # qseqid_syntenic_region = columns[0].split(";")[-1].split("::")[0]
                # sseqid_syntenic_region = columns[1].split("::")[0]
                #
                # if qseqid_syntenic_region == sseqid_syntenic_region:
                #    print(line.strip())
                qseqid_syntenic_region = columns[0].split(";")[-1].split("::")[0]
                qref = columns[0].split("::")[0]
                sseqid_syntenic_region = columns[1].split("::")[0]
                sseqid_info = columns[1].split("::")[1]

                if qseqid_syntenic_region == sseqid_syntenic_region:
                    passing_blast_lines.append(line.strip())

                    chrom, coord_info = sseqid_info.split(":")
                    start, end = map(int, coord_info.split("-"))

                    # Update the coordinates based on sstart and send
                    sstart, send = int(columns[8]), int(columns[9])

                    if sstart < send:
                        acr_start = start + sstart
                        acr_end = start + send
                    elif sstart > send:
                        acr_start = start + send
                        acr_end = start + sstart

                    final_name = "RefFrom;" + qref

                    print(final_name)

                    query_bed_line = f"{chrom}\t{acr_start}\t{acr_end}\t{final_name}\n"
                    passing_bed_lines.append(query_bed_line)

                    qref_region = columns[0].split("::")[1]
                    qref_chrom = qref_region.split(":")[0]
                    qref_start = qref_region.split(":")[1].split("-")[0]
                    qref_end = qref_region.split(":")[1].split("-")[1]
                    qref_bed_line = f"{qref_chrom}\t{qref_start}\t{qref_end}\t{qref}\n"
                    passing_ref_bed_lines.append(qref_bed_line)



            else:
                pass

    bed_file_output = output_file + ".blast_passing_regions.intersecting_regions.bed"
    bed_file_output_ref = output_file + ".blast_passing_regions.intersecting_regions.ref.bed"
    blast_passing_regions = output_file + ".blast_passing_regions.blast"

    remove_if_exists(bed_file_output)
    remove_if_exists(blast_passing_regions)

    with open(bed_file_output, 'a') as out:
        for line in passing_bed_lines:
            out.write(line)

    with open(blast_passing_regions, 'a') as out:
        for line in passing_blast_lines:
            out.write(line + "\n")
    with open(bed_file_output_ref, 'a') as out:
        for line in passing_ref_bed_lines:
            out.write(line)


#def get_parser():
#    parser = argparse.ArgumentParser(description='Calculate regions with fold \
#            enrichment over X and outputs these regions.')
#    parser.add_argument("-blast", help="Path to the BLAST format6 file.",
#                        required=True, dest="blast")
#    parser.add_argument("-evalue", type=float, default=0.001,
#                        dest="eval", help="Evalue threshold for filtering. Defaults to 0.001.")
#    parser.add_argument('-o', "--output", help=" Output file. If not \
#    given will be print ", required=True, dest='o')
#    args = vars(parser.parse_args())
#    return parser


#if __name__ == "__main__":
#    args = get_parser().parse_args()

#    print("Working on %s" % args.blast)
#    process_blast_file(args.blast, args.eval, args.o)




