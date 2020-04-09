# plot range of linear regressions for linear solver study

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from averageTime import *
from LinearRegression import *


def Scatter(solAlg, nonLinSol, variation):

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


    # find range for every linear solver
    column_names_Dense = column_names[0:7]
    y_axis_interception_Dense = y_axis_interceptions[0:7]
    slope_Dense = slopes[0:7]
    column_names_GMRES = column_names[7:14]
    y_axis_interception_GMRES = y_axis_interceptions[7:14]
    slope_GMRES = slopes[7:14]
    column_names_BiCGStab = column_names[14:21]
    y_axis_interception_BiCGStab = y_axis_interceptions[14:21]
    slope_BiCGStab = slopes[14:21]
    column_names_TFQMR = column_names[21:28]
    y_axis_interception_TFQMR = y_axis_interceptions[21:28]
    slope_TFQMR = slopes[21:28]
    column_names_KLU = column_names[28:35]
    y_axis_interception_KLU = y_axis_interceptions[28:35]
    slope_KLU = slopes[28:35]

    ax = plt.axes()
    if variation == 1:
        # sort values by slope to get boundaries and x-data
        column_names_Dense = [x for _,x in sorted(zip(slope_Dense, column_names_Dense))]
        y_axis_interception_Dense = [x for _,x in sorted(zip(slope_Dense, y_axis_interception_Dense))]
        slope_Dense = list(sorted(slope_Dense))

        column_names_GMRES = [x for _, x in sorted(zip(slope_GMRES, column_names_GMRES))]
        y_axis_interception_GMRES = [x for _,x in sorted(zip(slope_GMRES, y_axis_interception_GMRES))]
        slope_GMRES = list(sorted(slope_GMRES))

        column_names_BiCGStab = [x for _, x in sorted(zip(slope_BiCGStab, column_names_BiCGStab))]
        y_axis_interception_BiCGStab = [x for _,x in sorted(zip(slope_BiCGStab, y_axis_interception_BiCGStab))]
        slope_BiCGStab = list(sorted(slope_BiCGStab))

        column_names_TFQMR = [x for _, x in sorted(zip(slope_TFQMR, column_names_TFQMR))]
        y_axis_interception_TFQMR = [x for _,x in sorted(zip(slope_TFQMR, y_axis_interception_TFQMR))]
        slope_TFQMR = list(sorted(slope_TFQMR))

        column_names_KLU = [x for _, x in sorted(zip(slope_KLU, column_names_KLU))]
        y_axis_interception_KLU = [x for _,x in sorted(zip(slope_KLU, y_axis_interception_KLU))]
        slope_KLU = list(sorted(slope_KLU))

        ax.set_title('Sorted by slope')

    elif variation == 2:
        # sort values by y-axis-interception to get boundaries and x-data
        column_names_Dense = [x for _,x in sorted(zip(y_axis_interception_Dense, column_names_Dense))]
        slope_Dense = [x for _,x in sorted(zip(y_axis_interception_Dense, slope_Dense))]
        y_axis_interception_Dense = list(sorted(y_axis_interception_Dense))

        column_names_GMRES = [x for _, x in sorted(zip(y_axis_interception_GMRES, column_names_GMRES))]
        slope_GMRES = [x for _,x in sorted(zip(y_axis_interception_GMRES, slope_GMRES))]
        y_axis_interception_GMRES = list(sorted(y_axis_interception_GMRES))

        column_names_BiCGStab = [x for _, x in sorted(zip(y_axis_interception_BiCGStab, column_names_BiCGStab))]
        slope_BiCGStab = [x for _,x in sorted(zip(y_axis_interception_BiCGStab, slope_BiCGStab))]
        y_axis_interception_BiCGStab = list(sorted(y_axis_interception_BiCGStab))

        column_names_TFQMR = [x for _, x in sorted(zip(y_axis_interception_TFQMR, column_names_TFQMR))]
        slope_TFQMR = [x for _,x in sorted(zip(y_axis_interception_TFQMR, slope_TFQMR))]
        y_axis_interception_TFQMR = list(sorted(y_axis_interception_TFQMR))

        column_names_KLU = [x for _, x in sorted(zip(y_axis_interception_KLU, column_names_KLU))]
        slope_KLU = [x for _,x in sorted(zip(y_axis_interception_KLU, slope_KLU))]
        y_axis_interception_KLU = list(sorted(y_axis_interception_KLU))

        ax.set_title('Sorted by y-axis-interception')

    elif variation == 3:
        # sort values by both independent values to get boundaries and x-data
        column_names_Dense = [x for _,x in sorted(zip(slope_Dense, column_names_Dense))]
        y_axis_interception_Dense = list(sorted(y_axis_interception_Dense))
        slope_Dense = list(sorted(slope_Dense))

        column_names_GMRES = [x for _, x in sorted(zip(slope_GMRES, column_names_GMRES))]
        y_axis_interception_GMRES = list(sorted(y_axis_interception_GMRES))
        slope_GMRES = list(sorted(slope_GMRES))

        column_names_BiCGStab = [x for _, x in sorted(zip(slope_BiCGStab, column_names_BiCGStab))]
        y_axis_interception_BiCGStab = list(sorted(y_axis_interception_BiCGStab))
        slope_BiCGStab = list(sorted(slope_BiCGStab))

        column_names_TFQMR = [x for _, x in sorted(zip(slope_TFQMR, column_names_TFQMR))]
        y_axis_interception_TFQMR = list(sorted(y_axis_interception_TFQMR))
        slope_TFQMR = list(sorted(slope_TFQMR))

        column_names_KLU = [x for _, x in sorted(zip(slope_KLU, column_names_KLU))]
        y_axis_interception_KLU = list(sorted(y_axis_interception_KLU))
        slope_KLU = list(sorted(slope_KLU))

        ax.set_title('Sorted by slope and y-axis-interception independently --- data of slope')

    elif variation == 4:
        # sort values by both independent values to get boundaries and x-data
        column_names_Dense = [x for _,x in sorted(zip(y_axis_interception_Dense, column_names_Dense))]
        y_axis_interception_Dense = list(sorted(y_axis_interception_Dense))
        slope_Dense = list(sorted(slope_Dense))

        column_names_GMRES = [x for _, x in sorted(zip(y_axis_interception_GMRES, column_names_GMRES))]
        y_axis_interception_GMRES = list(sorted(y_axis_interception_GMRES))
        slope_GMRES = list(sorted(slope_GMRES))

        column_names_BiCGStab = [x for _, x in sorted(zip(y_axis_interception_BiCGStab, column_names_BiCGStab))]
        y_axis_interception_BiCGStab = list(sorted(y_axis_interception_BiCGStab))
        slope_BiCGStab = list(sorted(slope_BiCGStab))

        column_names_TFQMR = [x for _, x in sorted(zip(y_axis_interception_TFQMR, column_names_TFQMR))]
        y_axis_interception_TFQMR = list(sorted(y_axis_interception_TFQMR))
        slope_TFQMR = list(sorted(slope_TFQMR))

        column_names_KLU = [x for _, x in sorted(zip(y_axis_interception_KLU, column_names_KLU))]
        y_axis_interception_KLU = list(sorted(y_axis_interception_KLU))
        slope_KLU = list(sorted(slope_KLU))

        #ax.set_title('Sorted by slope and y-axis-interception independently --- data of y-axis-interception')

    # length of the last file
    file_length = len(next_tsv['id'])

    # get correct data for all five linear solvers and their linear regressions
    # plot a customized plot
    fontsize = 18
    labelsize = 12
    titlesize = 22

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

    colors = ['#d73027', '#fc8d59', '#fee090', '#91bfdb', '#4575b4']

    # num_x
    first_x_lowerbound = list(all_intern_columns[column_names.index(column_names_Dense[0])]['state_variables'])
    first_x_upperbound = list(all_intern_columns[column_names.index(column_names_Dense[6])]['state_variables'])
    second_x_lowerbound = list(all_intern_columns[column_names.index(column_names_GMRES[0]) + 7]['state_variables'])
    second_x_upperbound = list(all_intern_columns[column_names.index(column_names_GMRES[6]) + 7]['state_variables'])
    third_x_lowerbound = list(all_intern_columns[column_names.index(column_names_BiCGStab[0]) + 14]['state_variables'])
    third_x_upperbound = list(all_intern_columns[column_names.index(column_names_BiCGStab[6]) + 14]['state_variables'])
    fourth_x_lowerbound = list(all_intern_columns[column_names.index(column_names_TFQMR[0]) + 21]['state_variables'])
    fourth_x_upperbound = list(all_intern_columns[column_names.index(column_names_TFQMR[6]) + 21]['state_variables'])
    fifth_x_lowerbound = list(all_intern_columns[column_names.index(column_names_KLU[0]) + 28]['state_variables'])
    fifth_x_upperbound = list(all_intern_columns[column_names.index(column_names_KLU[6]) + 28]['state_variables'])

    # data
    first_data_lowerbound = [l[0] for l in [10 ** k for k in [y_axis_interception_Dense[0] + j for j in [slope_Dense[0] * i for i in first_x_lowerbound]]]]
    first_data_upperbound = [l[0] for l in [10 ** k for k in [y_axis_interception_Dense[6] + j for j in [slope_Dense[6] * i for i in first_x_upperbound]]]]
    second_data_lowerbound = [l[0] for l in [10 ** k for k in [y_axis_interception_GMRES[0] + j for j in [slope_GMRES[0] * i for i in second_x_lowerbound]]]]
    second_data_upperbound = [l[0] for l in [10 ** k for k in [y_axis_interception_GMRES[6] + j for j in [slope_GMRES[6] * i for i in second_x_upperbound]]]]
    third_data_lowerbound = [l[0] for l in [10 ** k for k in [y_axis_interception_BiCGStab[0] + j for j in [slope_BiCGStab[0] * i for i in third_x_lowerbound]]]]
    third_data_upperbound = [l[0] for l in [10 ** k for k in [y_axis_interception_BiCGStab[6] + j for j in [slope_BiCGStab[6] * i for i in third_x_upperbound]]]]
    fourth_data_lowerbound = [l[0] for l in [10 ** k for k in [y_axis_interception_TFQMR[0] + j for j in [slope_TFQMR[0] * i for i in fourth_x_lowerbound]]]]
    fourth_data_upperbound = [l[0] for l in [10 ** k for k in [y_axis_interception_TFQMR[6] + j for j in [slope_TFQMR[6] * i for i in fourth_x_upperbound]]]]
    fifth_data_lowerbound = [l[0] for l in [10 ** k for k in [y_axis_interception_KLU[0] + j for j in [slope_KLU[0] * i for i in fifth_x_lowerbound]]]]
    fifth_data_upperbound = [l[0] for l in [10 ** k for k in [y_axis_interception_KLU[6] + j for j in [slope_KLU[6] * i for i in fifth_x_upperbound]]]]

    # 10**num_x
    exp_first_x_lowerbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_Dense[0])]['state_variables'])]
    exp_first_x_upperbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_Dense[6])]['state_variables'])]
    exp_second_x_lowerbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_GMRES[0]) + 7]['state_variables'])]
    exp_second_x_upperbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_GMRES[6]) + 7]['state_variables'])]
    exp_third_x_lowerbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_BiCGStab[0]) + 14]['state_variables'])]
    exp_third_x_upperbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_BiCGStab[6]) + 14]['state_variables'])]
    exp_fourth_x_lowerbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_TFQMR[0]) + 21]['state_variables'])]
    exp_fourth_x_upperbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_TFQMR[6]) + 21]['state_variables'])]
    exp_fifth_x_lowerbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_KLU[0]) + 28]['state_variables'])]
    exp_fifth_x_upperbound = [10 ** m for m in list(all_intern_columns[column_names.index(column_names_KLU[6]) + 28]['state_variables'])]

    # change .tsv-id form e.g. 1_06_10.tsv to 06_10
    #linSol4legend_1 = 'DENSE'
    #linSol4legend_2 = 'GMRES '
    #linSol4legend_3 = 'BICGSTAB'
    #linSol4legend_4 = 'TFQMR'
    #linSol5legend_5 = 'KLU'

    # scatter plot
    ax.set_xlim([0.8, 1500])
    ax.set_ylim([0.1, 100000])  # 50000
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(exp_first_x_lowerbound, first_data_lowerbound, c=colors[0])#, label=str(linSol4legend_1))
    ax.plot(exp_first_x_upperbound, first_data_upperbound, c=colors[0])
    ax.plot(exp_second_x_lowerbound, second_data_lowerbound, c=colors[1])#, label=str(linSol4legend_2))
    ax.plot(exp_second_x_upperbound, second_data_upperbound, c=colors[1])
    ax.plot(exp_third_x_lowerbound, third_data_lowerbound, c=colors[2])#, label=str(linSol4legend_3))
    ax.plot(exp_third_x_upperbound, third_data_upperbound, c=colors[2])
    ax.plot(exp_fourth_x_lowerbound, fourth_data_lowerbound, c=colors[3])#, label=str(linSol4legend_4))
    ax.plot(exp_fourth_x_upperbound, fourth_data_upperbound, c=colors[3])
    ax.plot(exp_fifth_x_lowerbound, fifth_data_lowerbound, c=colors[4])#, label=str(linSol5legend_5))
    ax.plot(exp_fifth_x_upperbound, fifth_data_upperbound, c=colors[4])
    plt.tick_params(labelsize=labelsize)
    ax.set_xlabel('Number of state variables', fontsize=fontsize)
    ax.set_ylabel('Simulation time [ms]', fontsize=fontsize)
    # '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'

    # fill area under curves
    ax.fill([np.nanmin(exp_first_x_lowerbound), np.nanmax(exp_first_x_lowerbound), np.nanmax(exp_first_x_upperbound), np.nanmin(exp_first_x_upperbound)],
            [np.nanmin(first_data_lowerbound), np.nanmax(first_data_lowerbound), np.nanmax(first_data_upperbound), np.nanmin(first_data_upperbound)],
            c=colors[0], alpha=0.5, label='DENSE')
    ax.fill([np.nanmin(exp_second_x_lowerbound), np.nanmax(exp_second_x_lowerbound), np.nanmax(exp_second_x_upperbound), np.nanmin(exp_second_x_upperbound)],
            [np.nanmin(second_data_lowerbound), np.nanmax(second_data_lowerbound), np.nanmax(second_data_upperbound), np.nanmin(second_data_upperbound)],
            c=colors[1], alpha=0.5, label='GMRES')
    ax.fill([np.nanmin(exp_third_x_lowerbound), np.nanmax(exp_third_x_lowerbound), np.nanmax(exp_third_x_upperbound), np.nanmin(exp_third_x_upperbound)],
            [np.nanmin(third_data_lowerbound), np.nanmax(third_data_lowerbound), np.nanmax(third_data_upperbound), np.nanmin(third_data_upperbound)],
            c=colors[2], alpha=0.5, label='BICGSTAB')
    ax.fill([np.nanmin(exp_fourth_x_lowerbound), np.nanmax(exp_fourth_x_lowerbound), np.nanmax(exp_fourth_x_upperbound), np.nanmin(exp_fourth_x_upperbound)],
            [np.nanmin(fourth_data_lowerbound), np.nanmax(fourth_data_lowerbound), np.nanmax(fourth_data_upperbound), np.nanmin(fourth_data_upperbound)],
            c=colors[3], alpha=0.5, label='TFQMR')
    ax.fill([np.nanmin(exp_fifth_x_lowerbound), np.nanmax(exp_fifth_x_lowerbound), np.nanmax(exp_fifth_x_upperbound), np.nanmin(exp_fifth_x_upperbound)],
            [np.nanmin(fifth_data_lowerbound), np.nanmax(fifth_data_lowerbound), np.nanmax(fifth_data_upperbound), np.nanmin(fifth_data_upperbound)],
            c=colors[4], alpha=0.5, label='KLU')
    #ax.fill(np.append(exp_first_x_lowerbound, exp_first_x_upperbound[::-1]), np.append(first_data_lowerbound, first_data_upperbound[::-1]), c='#66c2a5', alpha=0.5)
    #ax.fill(np.append(exp_second_x_lowerbound, exp_second_x_upperbound[::-1]), np.append(second_data_lowerbound, second_data_upperbound[::-1]), c='#fc8d62', alpha=0.5)
    #ax.fill(np.append(exp_third_x_lowerbound, exp_third_x_upperbound[::-1]), np.append(third_data_lowerbound, third_data_upperbound[::-1]), c='#8da0cb', alpha=0.5)
    #ax.fill(np.append(exp_fourth_x_lowerbound, exp_fourth_x_upperbound[::-1]), np.append(fourth_data_lowerbound, fourth_data_upperbound[::-1]), c='#e78ac3', alpha=0.5)
    #ax.fill(np.append(exp_fifth_x_lowerbound, exp_fifth_x_upperbound[::-1]), np.append(fifth_data_lowerbound, fifth_data_upperbound[::-1]), c='#a6d854', alpha=0.5)


    # make top and right boxlines invisible
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # plot a legend
    ax.legend(loc=2, fontsize=labelsize, frameon=False)

    # better layout
    plt.tight_layout()

    # change plotting size
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)

    # save figure
    # plt.savefig('../paper_SolverSettings/Figures/Study_3/13012020/LinSol_' + solAlg + '_' + nonLinSol + '_Scatter.pdf')

    # show figure
    plt.show()


# call function
Scatter('2', '2', 4)