#!/usr/bin/env python

##ScAClass aims to predict the cell type identities

import os
import argparse
import sys
import subprocess
import re
import glob


##updating 110424 we will set the tissue and organ to be the reference
##updating 102224 we will split the training into different part to run the data to reduce the time by setting a low mem option
##updating 090723 add a step to update the clusterings

##preinstalled tools
##R: RcppML, harmony, Seurat, magrittr, matrixStats, Matrix
##python: smote, keras, sklearn, pandas


def get_parsed_args():

    parser = argparse.ArgumentParser(description="ScAClass pipeline")

    ##require files
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files."
                                                                    "Default: ./ ")

    parser.add_argument("-script_dir", dest='required_script_dir',
                        help='Users must provide the required script dir provided by GitHub.')

    parser.add_argument("-train_mtx", dest='train_acc_matrix', help="Provide the training matrix file downloaded from ScAClass.")

    parser.add_argument("-train_meta", dest='train_meta_file', help='Provide the meta file downloaded from ScAClass.')

    parser.add_argument("-species", dest = 'species_name', help='Provide a species name. Currently, the ScAClass has rice, maize, etc.')

    parser.add_argument("-organ", dest = 'organ_name', help = 'Provide a organ name.')

    parser.add_argument("-unknown_mtx", dest='unknown_acc_matrix', help = "Provide a matrix with unknown cell identity.")


    ##Optional parameters
    parser.add_argument("-open_low_mem",dest = 'open_low_mem', help = 'Users could use the low mem to perform the classification by selecting the top variant features of the training dataset.')

    parser.add_argument("-top_ft_num", dest = 'top_variant_ft_num', help = 'Once the low mem open, users could set a top variant ft number.')


    parser.add_argument("-clustering_fl", dest='clustering_fl_path', help = 'Provide a clustering file to build a clustering based annotation.')

    parser.add_argument("-running_steps", dest='running_steps', help = 'Provide running steps to check whether we will run the specific steps.'
                                                                       'Format is -running_steps 1,2,3,4'
                                                                       'Default: All')

    parser.add_argument("-PC_str" ,dest='PC_string', help = 'Provide a PC number string used for the embeddings, such as 30,50,70.'
                                                            'If users provide one more PC number, all the prediction under each PC number will be analyzed.'
                                                            'Default: 50.')

    parser.add_argument("-Rd_type", dest='Reduction_type', help = 'Provide a reduction type. Users can choose SVD or NMF.'
                                                                  'Default: SVD.')

    parser.add_argument("-ML_method", dest='ML_method', help = 'Provide a machine learning method. ScAClass provides RF (Random Forest) and SVM (Support Vector Machine).'
                                                               'Default: SVM.')

    parser.add_argument("-SVM_kernel", dest="SVM_kernel", help = 'SVM method provides ‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’ in sklearn.'
                                                                 'Default: rbf.')

    parser.add_argument("-open_smote", dest="open_smote", help = 'Smote will impute data enabling data balance manner.'
                                                                 'Default: yes')



    ##parse of parameters
    args = parser.parse_args()
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    ##for the input files
    if args.required_script_dir is None:
        print('Cannot find required script dir, please provide the dir in \'-script_dir\' !')
        return
    else:
        input_required_scripts_dir = args.required_script_dir


    if args.train_acc_matrix is None:
        print('Cannot find training matrix file, please provide it')
        return
    else:
        try:
            file = open(args.train_acc_matrix, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the training matrix file!')
            return

    if args.unknown_acc_matrix is None:
        print('Cannot find unknown matrix file, please provide it')
        return
    else:
        try:
            file = open(args.unknown_acc_matrix, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the unknown matrix file!')
            return

    if args.train_meta_file is None:
        print('Cannot find a train meta file, please provide it')
        return
    else:
        try:
            file = open(args.train_meta_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the train meta file!')
            return

    if args.clustering_fl_path is not None:
        try:
            file = open(args.clustering_fl_path, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the clustering file!')
            return

    if args.species_name is None:
        print('Cannot find species name, please provide it')
        return

    if args.organ_name is None:
        print('Cannot find organ name, please provide it')
        return


    ###########################################
    ##create the working and output directories
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    ##obtain the path of utils
    #run_script_path = __file__
    #if '/' in run_script_path:
    #    mt = re.match('(.+)/.+',run_script_path)
    #   run_script_dir = mt.group(1)
    #    utils_dir = run_script_dir + '/utils'

        ##running script dir needs to be change
        #utils_dir = '../input_other_required_scripts_dir/utils_ScAClass'

    #else:
    #    utils_dir = './utils'

    utils_dir = input_required_scripts_dir + '/utils_ScAClass'

    ##for the parameters
    if args.PC_string is None:
        PC_str = '50'
    else:
        PC_str = args.PC_string

    if args.Reduction_type is None:
        reduction_major_type = 'SVD'
    else:
        reduction_major_type = args.Reduction_type

    if args.ML_method is None:
        Machine_learning_method = 'svm'
    else:
        Machine_learning_method = args.ML_method

    if args.SVM_kernel is None:
        SVM_kernel_method = 'rbf'
    else:
        SVM_kernel_method = args.SVM_kernel

    if args.open_smote is None:
        decide_use_smote = 'yes'
    else:
        if args.open_smote == 'yes':
            decide_use_smote = 'yes'
        else:
            decide_use_smote = 'no'
            print('Users close the smote, please use \'-open_smote\' yes to open the smote')

    if args.running_steps is None:

        if args.clustering_fl_path is not None:
            final_running_steps = '1,2,3,4'
            final_running_steps_list = final_running_steps.split(',')
        else:
            final_running_steps = '1,2,3'
            final_running_steps_list = final_running_steps.split(',')
    else:
        final_running_steps = args.running_steps
        final_running_steps_list = final_running_steps.split(',')

    ##updating 102224
    ##default is not using the low mem
    if args.open_low_mem is None:
        decide_use_low_mem = 'no'
    else:
        if args.open_low_mem == 'yes':
            decide_use_low_mem = 'yes'
        else:
            decide_use_low_mem = 'no'

    if args.top_variant_ft_num is None:
        top_var_ft_num = '5000'
    else:
        top_var_ft_num = args.top_variant_ft_num


    ########
    ##Step01 prepare embedding files used for the training models
    step01_ID = '1'
    if step01_ID in final_running_steps_list:
        print('Step 1 prepare embedding files used for the training models')

        ##build the step01 directory
        Step01_prepare_embedding_dir = working_dir + '/Step01_prepare_embedding_dir'
        if not os.path.exists(Step01_prepare_embedding_dir):
            os.makedirs(Step01_prepare_embedding_dir)

        ipt_script = utils_dir + '/Step01_prepare_embeddings.py'

        ##updating 102224
        if decide_use_low_mem == 'no':
            print('Users decide not use the low mem')
            ipt_sub_script_1 = utils_dir + '/Step01_s1_feature_cell_embedings_updated.R'
        else:
            ipt_sub_script_1 = utils_dir + '/Step01_s1_feature_cell_embedings_updated_low_mem.R'

        input_training_matrix_rds_fl = args.train_acc_matrix
        input_unknown_matrix_csv_fl = args.unknown_acc_matrix
        input_feature_embedding_R_script_fl = ipt_sub_script_1

        ipt_train_meta_fl = args.train_meta_file

        cmd = 'python ' + ipt_script + \
              ' ' + input_training_matrix_rds_fl + \
              ' ' + input_unknown_matrix_csv_fl + \
              ' ' + input_feature_embedding_R_script_fl + \
              ' ' + ipt_train_meta_fl + \
              ' ' + PC_str + \
              ' ' + reduction_major_type + \
              ' ' + Step01_prepare_embedding_dir + \
              ' ' + top_var_ft_num
        print(cmd)
        subprocess.call(cmd,shell=True)

    ########
    ##Step02
    step02_ID = '2'
    if step02_ID in final_running_steps_list:
        print ('Step 2 training and make prediction')

        ##build the step02 directory
        Step02_training_prediction_dir = working_dir + '/Step02_training_prediction_dir'
        if not os.path.exists(Step02_training_prediction_dir):
            os.makedirs(Step02_training_prediction_dir)

        ipt_script = utils_dir + '/Step02_prediction.py'
        ipt_sub_script_1 = utils_dir + '/Step02_s1_prediction.py'

        ipt_train_meta_fl = args.train_meta_file
        cmd = 'python ' + ipt_script + \
              ' ' + ipt_sub_script_1 + \
              ' ' + ipt_train_meta_fl + \
              ' ' + working_dir + '/Step01_prepare_embedding_dir' + \
              ' ' + SVM_kernel_method + \
              ' ' + decide_use_smote + \
              ' ' + Machine_learning_method + \
              ' ' + Step02_training_prediction_dir + \
              ' ' + decide_use_low_mem
        print(cmd)
        subprocess.call(cmd,shell=True)

    ########
    ##step03
    step03_ID = '3'
    if step03_ID in final_running_steps_list:
        print ('Step 3 copy results to output dir')
        Step03_prepare_outputs_dir = working_dir + '/Step03_prepare_outputs_dir'
        if not os.path.exists(Step03_prepare_outputs_dir):
            os.makedirs(Step03_prepare_outputs_dir)

        store_final_line_list = []

        all_PC_dir_list = glob.glob(working_dir + '/Step02_training_prediction_dir' + '/*')
        for eachPC_dir in all_PC_dir_list:

            mt = re.match('.+/(.+)',eachPC_dir)
            PCdirnm = mt.group(1)

            mt = re.match('opt_PC(.+)_dir',PCdirnm)
            PCnum = mt.group(1)

            ##updating 102524
            if decide_use_low_mem == 'no':

                store_one_versus_all_dir = eachPC_dir + '/store_one_versus_all_dir'

                all_celltype_dir_list = glob.glob(store_one_versus_all_dir + '/*')

                for eachctdir in all_celltype_dir_list:

                    mt = re.match('.+/(.+)',eachctdir)
                    dirnm = mt.group(1)
                    mt = re.match('store_(.+)_dir',dirnm)
                    organct = mt.group(1)

                    opt_celltype_res_fl = eachctdir + '/' + '/opt_annotate_' + Machine_learning_method + '.txt'

                    with open (opt_celltype_res_fl,'r') as ipt:
                        for eachline in ipt:
                            eachline = eachline.strip('\n')
                            col = eachline.strip().split()
                            final_line = 'PC' + PCnum + '\t' + Machine_learning_method + '\t' + col[0] + '\t' + col[1]
                            store_final_line_list.append(final_line)

            ##udpating 102524
            else:
                all_range_dir_list = glob.glob(eachPC_dir + '/*')
                for eachrangdir in all_range_dir_list:

                    store_one_versus_all_dir = eachrangdir + '/store_one_versus_all_dir'

                    all_celltype_dir_list = glob.glob(store_one_versus_all_dir + '/*')

                    for eachctdir in all_celltype_dir_list:

                        mt = re.match('.+/(.+)', eachctdir)
                        dirnm = mt.group(1)
                        mt = re.match('store_(.+)_dir', dirnm)
                        organct = mt.group(1)

                        opt_celltype_res_fl = eachctdir + '/' + '/opt_annotate_' + Machine_learning_method + '.txt'

                        with open(opt_celltype_res_fl, 'r') as ipt:
                            for eachline in ipt:
                                eachline = eachline.strip('\n')
                                col = eachline.strip().split()
                                final_line = 'PC' + PCnum + '\t' + Machine_learning_method + '\t' + col[0] + '\t' + col[1]
                                store_final_line_list.append(final_line)

        with open (Step03_prepare_outputs_dir + '/opt_all_cells_predict_for_each_OneVersusAll_celltype.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        ##modify to each cell to have one line annotation
        store_PC_ml_cell_annotate_dic = {}
        with open (Step03_prepare_outputs_dir + '/opt_all_cells_predict_for_each_OneVersusAll_celltype.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                PC = col[0]
                method = col[1]
                cellnm = col[2]
                predict_celltype = col[3]

                if PC in store_PC_ml_cell_annotate_dic:
                    if method in store_PC_ml_cell_annotate_dic[PC]:
                        if cellnm in store_PC_ml_cell_annotate_dic[PC][method]:
                            store_PC_ml_cell_annotate_dic[PC][method][cellnm][predict_celltype] = 1
                        else:
                            store_PC_ml_cell_annotate_dic[PC][method][cellnm] = {}
                            store_PC_ml_cell_annotate_dic[PC][method][cellnm][predict_celltype] = 1
                    else:
                        store_PC_ml_cell_annotate_dic[PC][method] = {}
                        store_PC_ml_cell_annotate_dic[PC][method][cellnm] = {}
                        store_PC_ml_cell_annotate_dic[PC][method][cellnm][predict_celltype] = 1
                else:
                    store_PC_ml_cell_annotate_dic[PC] = {}
                    store_PC_ml_cell_annotate_dic[PC][method] = {}
                    store_PC_ml_cell_annotate_dic[PC][method][cellnm] = {}
                    store_PC_ml_cell_annotate_dic[PC][method][cellnm][predict_celltype] = 1

        ##write final results
        store_final_line_list = []
        first_line = 'PC' + '\t' + 'Method' + '\t' + 'Cell' + '\t' + 'Annotation'
        store_final_line_list.append(first_line)
        for eachPC in store_PC_ml_cell_annotate_dic:
            for eachmethod in store_PC_ml_cell_annotate_dic[eachPC]:
                for eachcell in store_PC_ml_cell_annotate_dic[eachPC][eachmethod]:
                    predict_celltype_dic = store_PC_ml_cell_annotate_dic[eachPC][eachmethod][eachcell]
                    if len(list(predict_celltype_dic.keys())) == 1:
                        annot = list(predict_celltype_dic.keys())[0]
                        if annot == 'Others':
                            final_annot = 'Unknown'
                        else:
                            final_annot = annot

                    else:
                        annot_list = list(predict_celltype_dic.keys())
                        new_annot_list = []
                        for eachitem in annot_list:
                            if eachitem != 'Others':
                                new_annot_list.append(eachitem)

                        final_annot = ','.join(new_annot_list)

                    final_line = eachPC + '\t' + eachmethod + '\t' + eachcell + '\t' + final_annot
                    store_final_line_list.append(final_line)

        with open (output_dir + '/opt_final_cell_annotation.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')


    ########
    ##step04
    ##build another cell type annotation based on the clustering results
    step04_ID = '4'
    if step04_ID in final_running_steps_list:

        print ('Step 4 use the clustering file to update the cell identities')
        ipt_script = utils_dir + '/Step04_update_cell_identity_by_cluster.py'
        cmd = 'python ' + ipt_script + \
              ' ' + args.clustering_fl_path + \
              ' ' + output_dir + '/opt_final_cell_annotation.txt' + \
              ' ' + output_dir
        print(cmd)
        subprocess.call(cmd,shell=True)



if __name__ == "__main__":
    main()















