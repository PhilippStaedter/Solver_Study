# script 2 to compare state trajectories from JWS with trajectories of the simulation
# => compares state trajectories

# Attention:    boundary conditions are not being simulated by JWS!

from execute_loadModels import *
from JWS_changeValues import *
from colourDataFrame import *
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import libsbml
import libsedml
import time
import statistics
import pandas as pd
import os
import urllib.request
import json
import itertools
import time


# upper and lower boundaries for the absolute and relative errors
AbsError_1 = range(-20,10)
RelError_2 = range(-20,10)

# create folder for all results
if not os.path.exists('./json_files_all_results_BDF'):                                                                # BDF
    os.makedirs('./json_files_all_results_BDF')

# iterate over all error combinations
for iAbsError in range(0, len(AbsError_1)):
    for iRelError in range(0, len(RelError_2)):

        if iAbsError != iRelError:
            continue

        # set errors
        abs_error = float('1e' + str(AbsError_1[iAbsError]))                                 # tighter conditions give back 'False' most of the time
        rel_error = float('1e' + str(RelError_2[iRelError]))

        # int2str
        abs_str = '{:.0e}'.format(float(abs_error))
        rel_str = '{:.0e}'.format(float(rel_error))
        
        print(f"TOLERANCES: abs={abs_str} rel={rel_str}")

        # create folder for all .csv files of the results
        if not os.path.exists('./json_files_all_results_BDF/json_files_' + abs_str + '_' + rel_str):
            os.makedirs('./json_files_all_results_BDF/json_files_' + abs_str + '_' + rel_str)

        # set counter
        counter = 0

        # get all models
        # list_directory_amici = sorted(os.listdir('../sbml2amici/amici_models_newest_version_0.10.19'))
        list_directory_amici = sorted(os.listdir('./sedml_models/'))

        # measure time needed for all mdoels
        start_time = time.time()

        # iterate over all models again
        for iMod in range(0, len(list_directory_amici)):
            iModel = list_directory_amici[iMod]
            #iModel = 'vanheerden2014_fig4-user'
            list_files = sorted(os.listdir('./sedml_models/' + iModel + '/sbml_models'))

            for iFile in list_files:
                print(f"    {iModel} :: {iFile}")

                # iFile without .sbml extension
                iFile, extension = iFile.split('.', 1)

                # important paths
                old_json_save_path = './json_files_BDF/' + iModel + '/' + iFile
                new_json_save_path = './json_files_all_results_BDF/json_files_' + abs_str + '_' + rel_str + '/' + iModel + '/' + iFile                    # BDF
                sedml_path = './sedml_models/' + iModel + '/' + iModel + '.sedml'
                sbml_path = './sedml_models/' + iModel + '/sbml_models/' + iFile + '.sbml'
                BioModels_path = './BioModelsDatabase_models'


                if os.path.exists(BioModels_path + '/' + iModel):
                    print('Model is still not part of JWS-database!')                                                                       # error 1
                else:

                    if not os.path.exists(old_json_save_path):
                        print('Model ' + iModel + '_' + iFile + ' crashed some other way!')                                                 # error 2
                    else:

                        # create folder
                        if not os.path.exists('./json_files_all_results_BDF/json_files_' + abs_str + '_' + rel_str + '/' + iModel + '/' + iFile):
                            os.makedirs('./json_files_all_results_BDF/json_files_' + abs_str + '_' + rel_str + '/' + iModel + '/' + iFile)

                        # open jws_simulation .csv file
                        tsv_file = pd.read_csv(old_json_save_path + '/' + iFile + '_JWS_simulation.csv', sep='\t')

                        # open model_simulation .csv file
                        df_state_trajectory = pd.read_csv(old_json_save_path + '/' + iFile + '_model_simulation.csv', sep='\t')

                        # columns names of .tsv file
                        column_names = list(tsv_file.columns)
                        column_names.remove('time')
                        del tsv_file['time']


                        # comparison
                        amount_col = len(column_names)
                        first_col = column_names[0]
                        amount_row = len(df_state_trajectory[first_col])
                        df_single_error = pd.DataFrame(columns=column_names, data=np.zeros((amount_row, amount_col)))
                        df_trajectory_error = pd.DataFrame(columns=column_names, data=np.zeros((1, amount_col)))
                        df_whole_error = pd.DataFrame(columns=['trajectories_match'], data=np.zeros((1, 1)))

                        # single error
                        for iCol in column_names:
                            for iRow in range(0, amount_row):
                                rel_tol = abs((df_state_trajectory.at[iRow, iCol] - tsv_file.at[iRow,iCol]) / df_state_trajectory.at[iRow, iCol])
                                abs_tol = abs(df_state_trajectory.at[iRow, iCol] - tsv_file.at[iRow,iCol])
                                if rel_tol <= rel_error or abs_tol <= abs_error:
                                    df_single_error.at[iRow, iCol] = 1
                                else:
                                    df_single_error.at[iRow, iCol] = 0

                        # trajectory error
                        for iCol in column_names:
                            if sum(df_single_error[iCol]) == amount_row:
                                df_trajectory_error.at[0, iCol] = 1
                            else:
                                df_trajectory_error.at[0, iCol] = 0

                        # whole error
                        error_list = []
                        for iCol in column_names:
                            error_list.append(df_trajectory_error.at[0, iCol])
                        if sum(error_list) == amount_col:
                            df_whole_error.at[0, 'trajectories_match'] = 1
                        else:
                            df_whole_error.at[0, 'trajectories_match'] = 0

                        # adjust counter
                        if df_whole_error.at[0, 'trajectories_match'] == 1:
                            print('matching state trajectory!')
                            counter = counter + 1

                        # save outcome
                        df_single_error.to_csv(path_or_buf=new_json_save_path + '/single_error.csv', sep='\t', index=False)
                        df_trajectory_error.to_csv(path_or_buf=new_json_save_path + '/trajectory_error.csv', sep='\t', index=False)
                        df_whole_error.to_csv(path_or_buf=new_json_save_path + '/whole_error.csv', sep='\t', index=False)

        # print number of all models with correct state trajectories
        print('Amount of models with correct state trajectories: ' + str(counter))

        print('time needed: ' + str(time.time() - start_time))
