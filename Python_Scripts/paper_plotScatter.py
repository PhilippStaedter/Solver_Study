# script to plot scatter plots for all 140 different combinations - study 3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from averageTime import *


def Scatter(solAlg, nonLinSol):

    # important paths
    base_path = '../paper_SolverSettings/WholeStudy'

    # list of all 35 data frames for better indexing in the future
    all_intern_columns = [pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]), pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[]),
                          pd.DataFrame(columns=[])]

    column_names = []

    # choose only the correct files
    all_files = sorted(os.listdir(base_path))
    correct_files = []
    for iFile in range(0, len(all_files)):
        if all_files[iFile].split('_')[0] == solAlg and all_files[iFile].split('_')[1] == nonLinSol:
            correct_files.append(all_files[iFile])

    # open all .tsv linear solver files + save right column in data frame
    for iCorrectFile in range(0, len(correct_files)):  # each .tsv file
        next_tsv = pd.read_csv(base_path + '/' + correct_files[iCorrectFile], sep='\t')

        # change .tsv-id form e.g. 1_06_10.tsv to 06_10
        new_name = correct_files[iCorrectFile].split('.')[0].split('_')[3] + '_' + correct_files[iCorrectFile].split('.')[0].split('_')[4]

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
        all_intern_columns[iCorrectFile]['state_variables'] = pd.Series(num_x)
        all_intern_columns[iCorrectFile][str(new_name)] = pd.Series(next_time_value)

    # length of the last file
    file_length = len(next_tsv['id'])

    # get correct data for all five linear solvers in one of the seven figures
    # plot a customized scatter plot
    fontsize = 22 - 12 + 2
    labelsize = 10 + 2
    titlesize = 30 - 8

    rotation = 90
    left = 0.07
    bottom = 0.75
    width = 0.4
    height = 0.18
    row_factor = 0.5
    column_factor = 0.22
    rotation_factor = 90
    alpha = 0.7

    linestyle = (0, (2, 5, 2, 5))
    linewidth = 0.1

    for iCounter in range(0, int(len(correct_files)/5)):

        # for ylim control
        if sorted(all_intern_columns[iCounter][column_names[iCounter]])[0] < 0.1 or sorted(all_intern_columns[iCounter + 7][column_names[iCounter + 7]])[0] < 0.1 or \
                sorted(all_intern_columns[iCounter + 14][column_names[iCounter + 14]])[0] < 0.1 or sorted(all_intern_columns[iCounter + 21][column_names[iCounter + 21]])[0] < 0.1 or \
                sorted(all_intern_columns[iCounter + 28][column_names[iCounter + 28]])[0] < 0.1:
            print('Need smaller ylim')

        # first plot
        if iCounter == 0:
            ax1 = plt.axes([left, bottom, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-6}$ + Rel. tol. = ' + r'$10^{-8}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            ax1.text(-0.1, 0.98, 'Simulation time [ms]', fontsize=fontsize, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 1:
            ax1 = plt.axes([left + iCounter * row_factor, bottom, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-8}$ + Rel. tol. = ' + r'$10^{-6}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            ax1.tick_params(labelleft=False)

        elif iCounter == 2:
            ax1 = plt.axes([left, bottom - (iCounter-1) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-8}$ + Rel. tol. = ' + r'$10^{-16}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            #ax1.tick_params(labelleft=False)
            ax1.text(-0.1, 0.98, 'Simulation time [ms]', fontsize=fontsize, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 3:
            ax1 = plt.axes([left + (iCounter-2) * row_factor, bottom - (iCounter-2) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-10}$ + Rel. tol. = ' + r'$10^{-12}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelleft=False)
            ax1.tick_params(labelbottom=False)
            #ax1.text(0.15, -0.3, 'Number of state variables', fontsize=fontsize, fontweight='bold', transform=ax1.transAxes)

        elif iCounter == 4:
            ax1 = plt.axes([left, bottom - (iCounter-2) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-12}$ + Rel. tol. = ' + r'$10^{-10}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.text(0.15, -0.3, 'Number of state variables', fontsize=fontsize, fontweight='bold', transform=ax1.transAxes)
            ax1.tick_params(labelbottom=False)
            ax1.text(-0.1, 0.98, 'Simulation time [ms]', fontsize=fontsize, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 5:
            ax1 = plt.axes([left + (iCounter-4) * row_factor, bottom - (iCounter-3) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-14}$ + Rel. tol. = ' + r'$10^{-14}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelleft=False)
            ax1.text(0.3, -0.3, 'Number of state variables', fontsize=fontsize, transform=ax1.transAxes)

        elif iCounter == 6:
            ax1 = plt.axes([left, bottom - (iCounter-3) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-16}$ + Rel. tol. = ' + r'$10^{-16}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.tick_params(labelbottom=False)ax1.set_ylim([0.1, 50000])
            ax1.text(0.3, -0.35, 'Number of state variables', fontsize=fontsize, transform=ax1.transAxes)
            ax1.text(-0.1, 0.98, 'Simulation time [ms]', fontsize=fontsize, transform=ax1.transAxes, rotation=rotation_factor)


        # apply formula:   iCounter --> iCounter + k*7
        # num_x
        first_x = all_intern_columns[iCounter]['state_variables']
        second_x = all_intern_columns[iCounter + 7]['state_variables']
        third_x = all_intern_columns[iCounter + 14]['state_variables']
        fourth_x = all_intern_columns[iCounter + 21]['state_variables']
        fifth_x = all_intern_columns[iCounter + 28]['state_variables']

        # data
        first_data = all_intern_columns[iCounter][column_names[iCounter]]
        second_data = all_intern_columns[iCounter + 7][column_names[iCounter + 7]]
        third_data = all_intern_columns[iCounter + 14][column_names[iCounter + 14]]
        fourth_data = all_intern_columns[iCounter + 21][column_names[iCounter + 21]]
        fifth_data = all_intern_columns[iCounter + 28][column_names[iCounter + 28]]

        # change .tsv-id form e.g. 1_06_10.tsv to 06_10
        linSol4legend_1 = 'Dense'
        linSol4legend_2 = 'SPGMR'
        linSol4legend_3 = 'SPBCG'
        linSol4legend_4 = 'SPTFQMR'
        linSol5legend_5 = 'KLU'

        # scatter plot
        ax1.set_xlim([0.8, 1500])
        ax1.set_ylim([0.1, 100000]) # 50000
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.scatter(first_x, first_data, alpha=alpha, c='#66c2a5', edgecolors='none', s=30, label=str(linSol4legend_1))# + ': ' + str(round(len(first_data)*100/file_length,2)) + ' %')
        ax1.scatter(second_x, second_data, alpha=alpha, c='#fc8d62', edgecolors='none', s=30, label=str(linSol4legend_2))# + ': ' + str(round(len(second_data)*100/file_length,2)) + ' %')
        ax1.scatter(third_x, third_data, alpha=alpha, c='#8da0cb', edgecolors='none', s=30, label=str(linSol4legend_3))# + ': ' + str(round(len(third_data)*100/file_length,2)) + ' %')
        ax1.scatter(fourth_x, fourth_data, alpha=alpha, c='#e78ac3', edgecolors='none', s=30, label=str(linSol4legend_4))# + ': ' + str(round(len(fourth_data)*100/file_length,2)) + ' %')
        ax1.scatter(fifth_x, fifth_data, alpha=alpha, c='#a6d854', edgecolors='none', s=30, label=str(linSol5legend_5))# + ': ' + str(round(len(fifth_data)*100/file_length,2)) + ' %')
        plt.tick_params(labelsize=labelsize)
        #'#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'

        # make top and right boxlines invisible
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # plot a legend
        if iCounter == 0:
            ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize, frameon=False)
        #if iCounter == 0:
        #    leg1 = ax1.legend(loc=4, frameon=False)
        '''
        leg2 = ax1.legend([str(round(len(all_intern_columns[iCounter]['state_variables'])*100/166,2)) + ' %',
                           str(round(len(all_intern_columns[iCounter + 1]['state_variables'])*100/166,2)) + ' %',
                           str(round(len(all_intern_columns[iCounter + 2]['state_variables'])*100/166,2)) + ' %',
                           str(round(len(all_intern_columns[iCounter + 3]['state_variables'])*100/166,2)) + ' %'], loc=4)
        ax1.add_artist(leg1)
        '''

    # plot a scatter plot in the right bottom corner
    alpha = 1
    rotation_factor = -30
    ax2 = plt.axes([left + row_factor, bottom - 3 * column_factor, width, height])
    for iCounter in range(0, int(len(correct_files) / 5)):
        first_x = 0.8 + iCounter
        second_x = 0.9 + iCounter
        third_x = 1 + iCounter
        fourth_x = 1.1 + iCounter
        fifth_x = 1.2 + iCounter
        first_data = all_intern_columns[iCounter][column_names[iCounter]]
        second_data = all_intern_columns[iCounter + 7][column_names[iCounter + 7]]
        third_data = all_intern_columns[iCounter + 14][column_names[iCounter + 14]]
        fourth_data = all_intern_columns[iCounter + 21][column_names[iCounter + 21]]
        fifth_data = all_intern_columns[iCounter + 28][column_names[iCounter + 28]]
        ax2.scatter(first_x, round(len(first_data)/file_length,2), alpha=alpha, c='#66c2a5', edgecolors='none', s=30,)
        ax2.scatter(second_x, round(len(second_data)/file_length,2), alpha=alpha, c='#fc8d62', edgecolors='none', s=30)
        ax2.scatter(third_x, round(len(third_data)/file_length,2), alpha=alpha, c='#8da0cb', edgecolors='none', s=30)
        ax2.scatter(fourth_x, round(len(fourth_data)/file_length,2), alpha=alpha, c='#e78ac3', edgecolors='none', s=30)
        ax2.scatter(fifth_x, round(len(fifth_data)/file_length,2), alpha=alpha, c='#a6d854', edgecolors='none', s=30)
        #ax2.axhline(round(len(fifth_data)/file_length,2) + iCounter * 0.01, fifth_x - 0.02, fifth_x + 0.02, c='#a6d854')
        plt.tick_params(labelsize=labelsize)
    ax2.set_ylim([0.7, 1.1])
    ax2.set_ylabel('Success rate [%]', fontsize=fontsize)
    ax2.set_yticklabels(['', '80', '90', '100'])

    # create major and minor ticklabels
    ax2.set_xticks([0, 1, 2, 3, 4, 5, 6, 7])
    ax2.set_xticks([0.001, 1.001, 2.001, 3.001, 4.001, 5.001, 6.001, 7.001], minor=True)
    ax2.set_xticklabels(['', r'$10^{-6}$', r'$10^{-8}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$'], fontsize=labelsize)
    ax2.set_xticklabels(['', r'$10^{-8}$', r'$10^{-6}$', r'$10^{-16}$', r'$10^{-12}$', r'$10^{-10}$', r'$10^{-14}$', r'$10^{-8}$'], minor=True, fontsize=labelsize)
    ax2.tick_params(axis='x', which='minor', pad=20)
    ax2.text(0.00001, -0.15, 'Abs. tol.: ', fontsize=fontsize, transform=ax2.transAxes)
    ax2.text(0.000025, -0.25, 'Rel. tol.: ', fontsize=fontsize, transform=ax2.transAxes)

    # make top and right boxlines invisible
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # set global labels
    #plt.text(0.1, 5.5, 'Simulation time distribution of models for different linear solver combinations', fontsize=titlesize, fontweight='bold', transform=ax1.transAxes)  # -60 , 350

    # better layout
    plt.tight_layout()

    # change plotting size
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)

    # save figure
    #plt.savefig('../paper_SolverSettings/Figures/Study_3/13012020/LinSol_' + solAlg + '_' + nonLinSol + '_Scatter.pdf')

    # show figure
    plt.show()


# call function
Scatter('2', '2')