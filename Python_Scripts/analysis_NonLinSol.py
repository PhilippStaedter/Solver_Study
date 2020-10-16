# Main Manuscript Plot --- Figure 2
# scrip to plot a Scatter and Histogram plot to investigate the difference between the non-linear solvers

import numpy as np
from averageTime import *
import math
import os
import pickle

data_file = "Data_NonLinSol_Main.pk"

# 1 is functional, 2 is Newton

def load_data():
    if os.path.exists(data_file):
        

    # check whether the folder 'Benchmarking_of_numerical_ODE_integration_methods/Data' exists
    if not os.path.exists('../../Benchmarking_of_numerical_ODE_integration_methods/Data/WholeStudy'):
        base_path = '../Data/WholeStudy'
    elif os.path.exists('../../Benchmarking_of_numerical_ODE_integration_methods/Data/WholeStudy'):
        base_path = '../../Benchmarking_of_numerical_ODE_integration_methods/Data/WholeStudy'

    # list of all data frames for nonLinSol == 1 for better indexing in the future
    all_intern_columns_1 = [pd.DataFrame(columns=[]) for _ in range(70)]

    # list of all data frames for nonLinSol == 2 for better indexing in the future
    all_intern_columns_2 = [pd.DataFrame(columns=[]) for _ in range(70)]

    all_intern_columns = [all_intern_columns_1, all_intern_columns_2]
    column_names = []

    # choose only the correct files
    all_files = sorted(os.listdir(base_path))
    correct_files_1 = []
    correct_files_2 = []
    for iFile in range(0, len(all_files)):
        if all_files[iFile].split('_')[0] == '1':
            correct_files_1.append(all_files[iFile])
        elif all_files[iFile].split('_')[0] == '2':
            correct_files_2.append(all_files[iFile])
    correct_files = [correct_files_1, correct_files_2]

    # open all .tsv linear solver files + save right column in data frame
    for iNonLinSol in range(0, len(correct_files)):
        for iCorrectFile in range(0, len(correct_files_1)):
            next_tsv = pd.read_csv(base_path + '/' + correct_files[iNonLinSol][iCorrectFile], sep='\t')

            # change .tsv-id form e.g. 1_06_10.tsv to 06_10
            new_name = correct_files[iNonLinSol][iCorrectFile].split('.')[0].split('_')[3] + '_' + \
                       correct_files[iNonLinSol][iCorrectFile].split('.')[0].split('_')[4]

            # reset after each iteration
            next_time_value = []
            num_x = []

            # open next file
            next_tsv = averaging(next_tsv)

            # get the correct values
            for iFile in range(0, len(next_tsv['id'])):  # each file
                if next_tsv['t_intern_ms'][iFile] != 0:
                    next_time_value.append(next_tsv['t_intern_ms'][iFile])
                    num_x.append(next_tsv['state_variables'][iFile])

            # append new column to existing data frame with correct values
            column_names.append(str(new_name))
            all_intern_columns[iNonLinSol][iCorrectFile]['state_variables'] = pd.Series(num_x)
            all_intern_columns[iNonLinSol][iCorrectFile][str(new_name)] = pd.Series(next_time_value)

    # length of the last file
    file_length = len(next_tsv['id'])

    # initialize y-data
    adams_data_1 = []
    bdf_data_1 = []
    adams_data_2 = []
    bdf_data_2 = []

    # Functional
    blank_space_counter = 0
    for iDensityPoint in range(0, int(len(correct_files_1)/2) + 4):
        if iDensityPoint in [7, 15, 23, 31]:
            adams_data_1.append(math.inf)
            bdf_data_1.append(math.inf)
            blank_space_counter += 1
        else:
            adams_data_1.append(1 - round(len(all_intern_columns_1[iDensityPoint - blank_space_counter][column_names[iDensityPoint - blank_space_counter]]) / file_length, 4))
            bdf_data_1.append(1 - round(len(all_intern_columns_2[iDensityPoint - blank_space_counter][column_names[iDensityPoint - blank_space_counter]]) / file_length, 4))

    # Newton-type
    blank_space_counter = 0
    for iDensityPoint in range(int(len(correct_files_1)/2), len(correct_files_1) + 4):
        if iDensityPoint in [int(len(correct_files_1)/2) + 7, int(len(correct_files_1)/2) + 15, int(len(correct_files_1)/2) + 23, int(len(correct_files_1)/2) + 31]:
            adams_data_2.append(math.inf)
            bdf_data_2.append(math.inf)
            blank_space_counter += 1
        else:
            adams_data_2.append(1 - round(len(all_intern_columns_1[iDensityPoint - blank_space_counter][column_names[iDensityPoint - blank_space_counter]]) / file_length, 4))
            bdf_data_2.append(1 - round(len(all_intern_columns_2[iDensityPoint - blank_space_counter][column_names[iDensityPoint - blank_space_counter]]) / file_length, 4))

    return adams_data_1, bdf_data_1, adams_data_2, bdf_data_2

# call functions
load_data()
