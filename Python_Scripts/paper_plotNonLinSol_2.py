# scrip to plot a bar plot to investigate the difference for the Multistep Method - study 5

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from averageTime import *
from matplotlib import ticker


def Multistep():

    # important paths
    base_path = '../paper_SolverSettings/WholeStudy'

    # list of all data frames for nonLinSol == 1 for better indexing in the future
    all_intern_columns_1 = [pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[])]

    # list of all data frames for nonLinSol == 2 for better indexing in the future
    all_intern_columns_2 = [pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[])]

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
        for iCorrectFile in range(0, len(correct_files_1)):  # each .tsv file
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

    # get correct data for the bar plot
    # plot a customized bar plot
    fontsize = 22 - 12 + 2
    labelsize = 10 + 4
    titlesize = 30 - 8

    rotation = 90
    left = 0.07
    bottom = 0.6
    width = 0.86
    height = 0.3
    row_factor = 0.45
    column_factor = 0.1
    rotation_factor = 70
    alpha = 1

    bar_width = 0.35

    # plot two bar plots
    figure = plt.figure()
    ax1 = figure.add_axes([left, bottom, width, height])
    ax2 = figure.add_axes([left, bottom - height - column_factor, width, height])
    index = np.arange(35)

    # just one for the legend
    adams_data_2 = all_intern_columns_1[int(len(correct_files_1)/2)][column_names[int(len(correct_files_1)/2)]]
    bdf_data_2 = all_intern_columns_2[int(len(correct_files_1)/2)][column_names[int(len(correct_files_1)/2)]]
    nonLinSol12 = ax2.bar(index[int(len(correct_files_1)/2) - int(len(correct_files_1) / 2)] - bar_width/2, round(len(adams_data_2) / file_length, 2), bar_width, alpha=alpha, color='orange', edgecolor='black', label='AM')
    nonLinSol22 = ax2.bar(index[int(len(correct_files_1)/2) - int(len(correct_files_1) / 2)] + bar_width/2, round(len(bdf_data_2) / file_length, 2), bar_width, alpha=alpha, color='blue', edgecolor='black', label='BDF')

    for iBarPlot in range(0, int(len(correct_files_1)/2)):
        adams_data_1 = all_intern_columns_1[iBarPlot][column_names[iBarPlot]]
        bdf_data_1 = all_intern_columns_2[iBarPlot][column_names[iBarPlot]]
        nonLinSol11 = ax1.bar(index[iBarPlot] - bar_width/2, round(len(adams_data_1)/file_length,2) , bar_width, alpha=alpha, edgecolor='black', color='orange')
        nonLinSol21 = ax1.bar(index[iBarPlot] + bar_width/2, round(len(bdf_data_1)/file_length,2), bar_width, alpha = alpha, edgecolor='black', color = 'blue')

    for iBarPlot in range(int(len(correct_files_1)/2) + 1, len(correct_files_1)):
        adams_data_2 = all_intern_columns_1[iBarPlot][column_names[iBarPlot]]
        bdf_data_2 = all_intern_columns_2[iBarPlot][column_names[iBarPlot]]
        nonLinSol12 = ax2.bar(index[iBarPlot - int(len(correct_files_1)/2)] - bar_width/2, round(len(adams_data_2)/file_length,2) , bar_width, alpha=alpha, edgecolor='black', color='orange')
        nonLinSol22 = ax2.bar(index[iBarPlot - int(len(correct_files_1)/2)] + bar_width/2, round(len(bdf_data_2)/file_length,2), bar_width, alpha = alpha, edgecolor='black', color = 'blue')

    ax1.set_title('Non-Linear solver: Functional', fontsize=titlesize)
    ax1.set_ylabel('Success rate [%]', fontsize=titlesize)
    ax2.set_title('Non-Linear solver: Newton-type', fontsize=titlesize)
    ax2.set_ylabel('Success rate [%]', fontsize=titlesize)
    ax1.set_xlim([-0.5, 34.5])
    ax2.set_xlim([-0.5, 34.5])
    ax1.set_ylim([0.7, 1])
    ax2.set_ylim([0.7, 1])
    ax1.set_xticklabels([])
    ax1.set_yticklabels(['70', '75', '80', '85', '90', '95', '100'], fontsize=fontsize)
    ax2.set_yticklabels(['70', '75', '80', '85', '90', '95', '100'], fontsize=fontsize)


    # create major and minor ticklabels
    upper_labels = ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'B', 'B', 'B', 'B', 'B', 'B',
                    'B', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'K', 'K', 'K', 'K', 'K', 'K', 'K']
    Abs_Rel_Tol = [r'$10^{-6}$' '\n' r'$10^{-8}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-8}$' '\n' r'$10^{-6}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-8}$' '\n' r'$10^{-16}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-10}$' '\n' r'$10^{-12}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-12}$' '\n'  r'$10^{-10}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-14}$' '\n' r'$10^{-14}$', ' \n ', ' \n ', ' \n ', ' \n ',
                   r'$10^{-16}$' '\n' r'$10^{-8}$', ' \n ', ' \n ', ' \n ', ' \n ']

    ax1.set_xticks(list(range(35)))
    ax2.set_xticks(list(range(35)))
    minor_list_1 = [x + 0.001 for x in list(range(35))]
    ax2.set_xticks(minor_list_1, minor=True)
    ax2.set_xticklabels(upper_labels, fontsize=fontsize)
    ax2.set_xticklabels(Abs_Rel_Tol, minor=True, fontsize=fontsize)
    ax2.tick_params(axis='x', which='major', pad=10)
    ax2.tick_params(axis='x', which='minor', pad=30)
    ax2.text(-0.05, -0.1, 'Lin. sol.: ', fontsize=fontsize, transform=ax2.transAxes)
    ax2.text(-0.05, -0.18, 'Abs. tol.: ', fontsize=fontsize, transform=ax2.transAxes)
    ax2.text(-0.05, -0.27, 'Rel. tol.: ', fontsize=fontsize, transform=ax2.transAxes)

    # create new empty invisible axis for legend
    #ax3 = figure.add_axes([0.15, 0.4, 0.02, 0.02])
    #ax3.plot(range(2), c='orange', label='AM')
    #ax3.plot(range(2), c='blue', label='BDF')
    ax2.legend(loc='upper left', bbox_to_anchor=(0, -0.3), fancybox=True, shadow=True, ncol=5, frameon=False, fontsize=titlesize)
    ax2.text(0.4, -0.5, 'D: Dense,  G: GMRES,  B: BCG,  T: TFQMR,  K: KLU', fontsize=titlesize, transform=ax2.transAxes)

    # make top and right boxlines invisible
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # set global labels
    #ax1.set_title('Adams-Moulton vs BDF - Success Rate', fontsize=titlesize, fontweight='bold', pad=30)
    #plt.text(0.3, 2.75, 'Adams-Moulton vs BDF - Success Rate', fontsize=titlesize, fontweight='bold', transform=ax2.transAxes)

    # better layout
    plt.tight_layout()

    # change plotting size
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)

    # save figure
    #plt.savefig('../paper_SolverSettings/Figures/Study_5/Success_Rate_SolAlg.pdf')

    # show figure
    plt.show()

# call functions
Multistep()