#!/usr/bin/env python

##this script is to train and predict the unknown cells
##updating 101225 we will try different methods of training

from sklearn.ensemble import RandomForestClassifier
from pandas import read_csv
from sklearn import svm
import pickle
import numpy as np
import sys
from imblearn.over_sampling import SMOTE

##updating 101225
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import GradientBoostingClassifier



##updating 123125
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline



input_training_dataset_fl = sys.argv[1]

input_indep_testing_dataset_fl = sys.argv[2]

input_training_meta_fl = sys.argv[3]

input_output_dir = sys.argv[4]

SVM_para = sys.argv[5]

open_smote = sys.argv[6]

target_method_string = sys.argv[7]

target_method_list = target_method_string.split(',')


##we will directly write out the outputs
def write_out_results (labels_pred,id_nm_dic,features_indep_test,testing_type,input_opt_dir):

    ##testing_type is SVMtest or SVMindetest RFtest RFindetest

    #################
    ##updating 042522
    pred_nm_list = []
    for eachpred_id in labels_pred[0]:
        # print(eachpred_id)
        pred_nm = id_nm_dic[str(eachpred_id)]
        pred_nm_list.append(pred_nm)

    ##make sure features_indep_test is the format of cell by features
    store_final_line_list = []
    list_count = -1
    for eachcellnm in features_indep_test.index:
        list_count += 1
        pred_celltypenm = pred_nm_list[list_count]
        final_line = eachcellnm + '\t' + pred_celltypenm
        store_final_line_list.append(final_line)

    with open(input_opt_dir + '/opt_annotate_' + testing_type + '.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def train_evaluation (input_train_file,input_meta_train_file,
                      input_test_indep_file,input_opt_dir,SVM_para,open_smote,target_method_list):

    ##run script
    ##load import file
    features_train = read_csv(input_train_file, index_col = 0).values.T
    meta_train = read_csv(input_meta_train_file,header = 0, index_col = 0)

    label_names = meta_train.loc[:,"cell_type"].unique()
    names_to_id = { key:val for val,key in enumerate(label_names)}
    labels_train = meta_train.loc[:,"cell_type"].replace(names_to_id).values

    ##updating 041522 decide whether we need the smote
    if open_smote == 'yes':
        print('we will use the smote to do the augmentation')
        sm = SMOTE(random_state=42)
        features_train, labels_train = sm.fit_resample(features_train, labels_train)

    ##initiate a function to store single cell id and its real name
    ##change the order of id and name
    id_nm_dic = {}
    for eachnm in names_to_id:
        id_nm_dic[str(names_to_id[eachnm])] = eachnm

    #########
    ##For svm
    ##PREIOUS is linear
    if 'svm' in target_method_list:

        ##updating 122125
        svm_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("svm", svm.SVC(probability=True))
        ])

        ##hyperparameter tunning
        param_grid = [
            # RBF kernel
            {
                "svm__kernel": ["rbf"],
                "svm__C": [0.1, 1, 10, 100],
                "svm__gamma": ["scale", 0.01, 0.1, 1]
            },

            # Linear kernel
            {
                "svm__kernel": ["linear"],
                "svm__C": [0.01, 0.1, 1, 10, 100]
            }
        ]

        grid = GridSearchCV(
            estimator=svm_pipe,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf_svm = grid.best_estimator_

        pkl_filename = input_opt_dir + "/svm_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf_svm, file)

        with open(input_opt_dir + "/svm_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))



        #clf_svm = svm.SVC(kernel=SVM_para,probability=True)
        #clf_svm.fit(features_train, labels_train)
        ##save svm model
        #pkl_filename = input_opt_dir + "/svm_model.pkl"
        #with open(pkl_filename, 'wb') as file:
        #    pickle.dump(clf_svm, file)

        ##for the validation data
        clfs = [clf_svm]

        ##for the independent data
        features_indep_test =  read_csv(input_test_indep_file, index_col = 0).values.T
        #meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        #labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
        features_indep_test_novalue = read_csv(input_test_indep_file, index_col=0).T

        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array([prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        ##close write compare file
        #writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest',only_predict_unknown_cell,'no')
        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_novalue, 'svm', input_opt_dir)


    ########
    ##For RF
    if 'rf' in target_method_list:


        rf = RandomForestClassifier(
            random_state=42,
            n_jobs=1
        )

        #param_grid = {
        #    "n_estimators": [200, 500, 1000],
        #    "max_depth": [None, 10, 20, 40],
        #    "min_samples_split": [2, 5, 10],
        #    "min_samples_leaf": [1, 2, 4],
        #    "max_features": ["sqrt", "log2"]
        #}

        param_grid = {
            "n_estimators": [300, 500],
            "max_depth": [10, 20, None],
            "min_samples_leaf": [1, 2, 5],
            "max_features": ["sqrt", "log2"]
        }


        grid = GridSearchCV(
            estimator=rf,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf_rf = grid.best_estimator_


        #clf_rf = RandomForestClassifier(n_estimators = 500, n_jobs = 42).fit(features_train, labels_train)
        ##save rf model
        #pkl_filename = input_opt_dir + "/rf_model.pkl"
        #with open(pkl_filename, 'wb') as file:
        #    pickle.dump(clf_rf, file)

        pkl_filename = input_opt_dir + "/rf_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf_rf, file)

        # （可选）保存最优参数，方便审稿
        with open(input_opt_dir + "/rf_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))


        ##for the validation data
        clfs = [clf_rf]

        ##for the independent data
        features_indep_test =  read_csv(input_test_indep_file, index_col = 0).values.T

        features_indep_test_notvalue = read_csv(input_test_indep_file, index_col=0).T
        #meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        #labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values

        ##updating 020221
        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array([prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_notvalue, 'rf', input_opt_dir)

    ##updating 101225
    #if 'others' in target_method_list:

    if 'lr' in target_method_list:

        lr_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("lr", LogisticRegression(
                max_iter=2000,
                multi_class="auto"
            ))
        ])

        param_grid = [
            {
                "lr__penalty": ["l2"],
                "lr__C": [0.01, 0.1, 1, 10, 100],
                "lr__solver": ["lbfgs", "saga"]
            },
            {
                "lr__penalty": ["l1"],
                "lr__C": [0.01, 0.1, 1, 10],
                "lr__solver": ["liblinear", "saga"]
            }
        ]

        grid = GridSearchCV(
            estimator=lr_pipe,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf = grid.best_estimator_

        pkl_filename = input_opt_dir + "/lr_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf, file)

        # （可选）保存最优参数，方便审稿
        with open(input_opt_dir + "/lr_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))



        #models = {
        #    "lr": LogisticRegression(max_iter=1000)
        #}

        #for name, model in models.items():
        #    clf = model.fit(features_train, labels_train)

        name = 'lr'

        #pkl_filename = input_opt_dir + "/" + name + "_model.pkl"
        #with open(pkl_filename, 'wb') as file:
        #    pickle.dump(clf, file)

        ##for the validation data
        clfs = [clf]

        ##for the independent data
        features_indep_test = read_csv(input_test_indep_file, index_col=0).values.T
        # meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        # labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
        features_indep_test_novalue = read_csv(input_test_indep_file, index_col=0).T

        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array(
                [prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        ##close write compare file
        # writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest',only_predict_unknown_cell,'no')
        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_novalue, name, input_opt_dir)

    if 'knn' in target_method_list:

        knn_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("knn", KNeighborsClassifier())
        ])

        param_grid = [
            {
                "knn__n_neighbors": [3, 5, 7, 11, 15, 21],
                "knn__weights": ["uniform", "distance"],
                "knn__metric": ["euclidean"]
            },
            {
                "knn__n_neighbors": [3, 5, 7, 11, 15, 21],
                "knn__weights": ["uniform", "distance"],
                "knn__metric": ["manhattan"]
            }
        ]

        grid = GridSearchCV(
            estimator=knn_pipe,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf = grid.best_estimator_

        pkl_filename = input_opt_dir + "/knn_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf, file)

        # （可选）保存最优参数，方便审稿
        with open(input_opt_dir + "/knn_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))


        #models = {
        #    "knn": KNeighborsClassifier(n_neighbors=5)
        #}

        #for name, model in models.items():
            #clf = model.fit(features_train, labels_train)

            #pkl_filename = input_opt_dir + "/" + name + "_model.pkl"
            #with open(pkl_filename, 'wb') as file:
            #    pickle.dump(clf, file)

            ##for the validation data
        clfs = [clf]

        ##for the independent data
        features_indep_test = read_csv(input_test_indep_file, index_col=0).values.T
        # meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        # labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
        features_indep_test_novalue = read_csv(input_test_indep_file, index_col=0).T

        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array(
                [prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        name = 'knn'

        ##close write compare file
        # writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest',only_predict_unknown_cell,'no')
        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_novalue, name, input_opt_dir)

    if 'nb' in target_method_list:

        nb = GaussianNB()

        param_grid = {
            "var_smoothing": [1e-12, 1e-10, 1e-9, 1e-8, 1e-7]
        }

        grid = GridSearchCV(
            estimator=nb,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf = grid.best_estimator_

        pkl_filename = input_opt_dir + "/nb_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf, file)

        # （可选）保存最优参数，方便写 Methods / Response
        with open(input_opt_dir + "/nb_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))

        #models = {
        #    "nb": GaussianNB()
        #}

        #for name, model in models.items():
            #clf = model.fit(features_train, labels_train)

            #pkl_filename = input_opt_dir + "/" + name + "_model.pkl"
            #with open(pkl_filename, 'wb') as file:
            #    pickle.dump(clf, file)

        ##for the validation data
        clfs = [clf]

        ##for the independent data
        features_indep_test = read_csv(input_test_indep_file, index_col=0).values.T
        # meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        # labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
        features_indep_test_novalue = read_csv(input_test_indep_file, index_col=0).T

        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array(
                [prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        name = 'nb'

        ##close write compare file
        # writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest',only_predict_unknown_cell,'no')
        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_novalue, name, input_opt_dir)

    if 'gbm' in target_method_list:

        gbm = GradientBoostingClassifier(random_state=42)

        param_grid = {
            "n_estimators": [200, 500],
            "learning_rate": [0.01, 0.05, 0.1],
            "max_depth": [3, 5],
            "subsample": [0.8, 1.0]
        }

        grid = GridSearchCV(
            estimator=gbm,
            param_grid=param_grid,
            cv=5,
            scoring="f1_macro",
            n_jobs=-1
        )

        grid.fit(features_train, labels_train)

        clf = grid.best_estimator_

        pkl_filename = input_opt_dir + "/gbm_model_best.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf, file)

        # （可选但强烈推荐）保存最优参数
        with open(input_opt_dir + "/gbm_best_params.txt", "w") as f:
            f.write(str(grid.best_params_))



        #models = {
        #    "gbm": GradientBoostingClassifier(n_estimators=500)
        #}

        #for name, model in models.items():
            #clf = model.fit(features_train, labels_train)

            #pkl_filename = input_opt_dir + "/" + name + "_model.pkl"
            #with open(pkl_filename, 'wb') as file:
            #    pickle.dump(clf, file)

        ##for the validation data
        clfs = [clf]

        ##for the independent data
        features_indep_test = read_csv(input_test_indep_file, index_col=0).values.T
        # meta_indep_test = read_csv(input_meta_test_indep_file,header = 0, index_col = 0)
        # labels_indep_test = meta_indep_test.loc[:,"cell_type"].replace(names_to_id).values
        features_indep_test_novalue = read_csv(input_test_indep_file, index_col=0).T

        prob_pred_all = ''
        labels_indep_pred = []
        prob_indep_pred = []
        for clf in clfs:
            prob_pred_all = clf.predict_proba(features_indep_test)
            labels_indep_pred.append(np.argsort(-prob_pred_all, axis=1)[:, 0])
            prob_indep_pred.append(np.array(
                [prob_pred_all[i, labels_indep_pred[-1][i]] for i in range(labels_indep_pred[-1].shape[0])]))

        name = 'gbm'

        ##close write compare file
        # writeout_compare_file(labels_indep_pred, labels_indep_test, meta_indep_test, id_nm_dic, input_opt_dir, 'SVMindetest',only_predict_unknown_cell,'no')
        write_out_results(labels_indep_pred, id_nm_dic, features_indep_test_novalue, name, input_opt_dir)




train_evaluation (input_training_dataset_fl,input_training_meta_fl,
                  input_indep_testing_dataset_fl,input_output_dir,SVM_para,open_smote,target_method_list)




