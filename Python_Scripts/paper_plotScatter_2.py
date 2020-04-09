# script to plot linear regressions for all scatter plots - study 3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from averageTime import *
from LinearRegression import *


def Scatter(solAlg, nonLinSol):

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

    # important paths
    base_path = '../paper_SolverSettings/WholeStudy'

    # choose only the correct files
    all_files = sorted(os.listdir(base_path))
    correct_files = []
    for iFile in range(0, len(all_files)):
        if all_files[iFile].split('_')[0] == solAlg and all_files[iFile].split('_')[1] == nonLinSol:
            correct_files.append(all_files[iFile])

    # open all .tsv linear solver files + save right column in data frame
    y_axis_interceptions = []
    slopes = []
    for iCorrectFile in range(0, len(correct_files)):  # each .tsv file
        next_tsv = pd.read_csv(base_path + '/' + correct_files[iCorrectFile], sep='\t')

        # change .tsv-id form e.g. 1_06_10.tsv to 06_10
        new_name = correct_files[iCorrectFile].split('.')[0].split('_')[3] + '_' + correct_files[iCorrectFile].split('.')[0].split('_')[4]

        # reset after each iteration
        next_time_value = []
        num_x = []

        # open next file
        next_tsv = averaging(next_tsv)

        # do a linear regression
        y_axis_interception, slope = linearRegression(next_tsv, 'state_variables', 't_intern_ms')
        y_axis_interceptions.append(y_axis_interception)
        slopes.append(slope)

        # get the correct values
        for iFile in range(0, len(next_tsv['id'])):  # each file
            if next_tsv['t_intern_ms'][iFile] != 0:
                next_time_value.append(np.log10(next_tsv['t_intern_ms'][iFile]))
                num_x.append(np.log10(next_tsv['state_variables'][iFile]))

        # append new column to existing data frame with correct log10 values
        column_names.append(str(new_name))
        all_intern_columns[iCorrectFile]['state_variables'] = pd.Series(num_x)
        all_intern_columns[iCorrectFile][str(new_name)] = pd.Series(next_time_value)


    # length of the last file
    file_length = len(next_tsv['id'])

    # get correct data for all five linear solvers in one of the seven figures
    # plot a customized plot
    fontsize = 7
    labelsize = 7
    titlesize = 22

    rotation = 90
    left = 0.07
    bottom = 0.75
    width = 0.44
    height = 0.18
    row_factor = 0.47
    column_factor = 0.22
    rotation_factor = 90
    alpha = 0.7

    linestyle = (0, (2, 5, 2, 5))
    linewidth = 0.1

    for iCounter in range(0, int(len(correct_files)/5)):

        # first plot
        if iCounter == 0:
            ax1 = plt.axes([left, bottom, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-6}$ + Rel. tol. = ' + r'$10^{-8}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            ax1.text(-0.12, 0.07, 'Simulation time [ms]', fontsize=fontsize - 1, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 1:
            ax1 = plt.axes([left + iCounter * row_factor, bottom, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-8}$ + Rel. tol. = ' + r'$10^{-6}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            ax1.tick_params(labelleft=False)

        elif iCounter == 2:
            ax1 = plt.axes([left, bottom - (iCounter - 1) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-8}$ + Rel. tol. = ' + r'$10^{-16}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            ax1.tick_params(labelbottom=False)
            #ax1.tick_params(labelleft=False)
            ax1.text(-0.12, 0.04, 'Simulation time [ms]', fontsize=fontsize - 1, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 3:
            ax1 = plt.axes([left + (iCounter - 2) * row_factor, bottom - (iCounter-2) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-10}$ + Rel. tol. = ' + r'$10^{-12}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelleft=False)
            ax1.tick_params(labelbottom=False)
            #ax1.text(0.15, -0.3, 'Number of state variables', fontsize=fontsize, fontweight='bold', transform=ax1.transAxes)

        elif iCounter == 4:
            ax1 = plt.axes([left, bottom - (iCounter - 2) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-12}$ + Rel. tol. = ' + r'$10^{-10}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.text(0.15, -0.3, 'Number of state variables', fontsize=fontsize, fontweight='bold', transform=ax1.transAxes)
            ax1.tick_params(labelbottom=False)
            ax1.text(-0.12, 0.01, 'Simulation time [ms]', fontsize=fontsize - 1, transform=ax1.transAxes, rotation=rotation_factor)

        elif iCounter == 5:
            ax1 = plt.axes([left + (iCounter - 4) * row_factor, bottom - (iCounter-3) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-14}$ + Rel. tol. = ' + r'$10^{-14}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_yaxis().set_visible(False)
            ax1.tick_params(labelleft=False)
            ax1.text(0.3, -0.32, 'Number of state variables', fontsize=fontsize - 1, transform=ax1.transAxes)

        elif iCounter == 6:
            ax1 = plt.axes([left, bottom - (iCounter - 3) * column_factor, width, height])
            ax1.text(0.3, 1.05, 'Abs. tol. = ' + r'$10^{-16}$ + Rel. tol. = ' + r'$10^{-16}$', fontsize=fontsize, transform=ax1.transAxes)
            #ax1.get_xaxis().set_visible(False)
            #ax1.tick_params(labelbottom=False)ax1.set_ylim([0.1, 50000])
            ax1.text(0.3, -0.38, 'Number of state variables', fontsize=fontsize - 1, transform=ax1.transAxes)
            ax1.text(-0.12, -0.02, 'Simulation time [ms]', fontsize=fontsize - 1, transform=ax1.transAxes, rotation=rotation_factor)


        # apply formula:   iCounter --> iCounter + k*7
        # num_x
        first_x = list(all_intern_columns[iCounter]['state_variables'])
        second_x = list(all_intern_columns[iCounter + 7]['state_variables'])
        third_x = list(all_intern_columns[iCounter + 14]['state_variables'])
        fourth_x = list(all_intern_columns[iCounter + 21]['state_variables'])
        fifth_x = list(all_intern_columns[iCounter + 28]['state_variables'])

        # data
        first_data = [l[0] for l in [10**k for k in [y_axis_interceptions[iCounter] + j for j in [slopes[iCounter]*i for i in first_x]]]]
        second_data = [l[0] for l in [10**k for k in [y_axis_interceptions[iCounter + 7] + j for j in [slopes[iCounter + 7]*i for i in second_x]]]]
        third_data = [l[0] for l in [10**k for k in [y_axis_interceptions[iCounter + 14] + j for j in [slopes[iCounter + 14]*i for i in third_x]]]]
        fourth_data = [l[0] for l in [10**k for k in [y_axis_interceptions[iCounter + 21] + j for j in [slopes[iCounter + 21]*i for i in fourth_x]]]]
        fifth_data = [l[0] for l in [10**k for k in [y_axis_interceptions[iCounter + 28] + j for j in [slopes[iCounter + 28]*i for i in fifth_x]]]]

        # 10**num_x
        exp_first_x = [10**m for m in list(all_intern_columns[iCounter]['state_variables'])]
        exp_second_x = [10**m for m in list(all_intern_columns[iCounter + 7]['state_variables'])]
        exp_third_x = [10**m for m in list(all_intern_columns[iCounter + 14]['state_variables'])]
        exp_fourth_x = [10**m for m in list(all_intern_columns[iCounter + 21]['state_variables'])]
        exp_fifth_x = [10**m for m in list(all_intern_columns[iCounter + 28]['state_variables'])]

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
        ax1.plot(exp_first_x, first_data, c='#66c2a5', label=str(linSol4legend_1) + '_slope: ' + str(slopes[iCounter]))
        ax1.plot(exp_second_x, second_data, c='#fc8d62', label=str(linSol4legend_2) + '_slope: ' + str(slopes[iCounter + 7]))
        ax1.plot(exp_third_x, third_data, c='#8da0cb', label=str(linSol4legend_3) + '_slope: ' + str(slopes[iCounter + 14]))
        ax1.plot(exp_fourth_x, fourth_data, c='#e78ac3', label=str(linSol4legend_4) + '_slope: ' + str(slopes[iCounter + 21]))
        ax1.plot(exp_fifth_x, fifth_data, c='#a6d854', label=str(linSol5legend_5) + '_slope: ' + str(slopes[iCounter + 28]))
        plt.tick_params(labelsize=labelsize)
        #'#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'

        # make top and right boxlines invisible
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # plot a legend
        ax1.legend(loc=2, fontsize=fontsize - 3, frameon=False)


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