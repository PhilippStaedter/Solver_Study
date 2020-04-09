# script to create box plot with percentiles and median - with success rates

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from averageTime import *
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


Multistep_Method = 'BDF'


# decision function
if Multistep_Method == 'Adams':
    prefix = '1_'
elif Multistep_Method == 'BDF':
    prefix = '2_'

# important paths
#tolerance_path = '../bachelor_thesis/Tolerance'
tolerance_path = '../paper_SolverSettings/Tolerances_1e4/' + Multistep_Method                                                           #

# main .tsv file to norm all other files
main_tsv = pd.read_csv(tolerance_path + '/' + prefix + '06_06.tsv', sep='\t')                                                       #

# get new .tsv file
main_tsv = averaging(main_tsv)

# set two axes objects
figure = plt.figure()
ax1 = figure.add_axes([0.1, 0.55, 0.8, 0.4])                 # ax = plt.axes()
ax2 = figure.add_axes([0.1, 0.15, 0.8, 0.35])

# get list for all data
xTickLabel = []
total_data = []

# open all .tsv tolerance files
tolerance_files_old = sorted(os.listdir(tolerance_path))
#del tolerance_files[0]
tolerance_files = []
for iTolFile in range(0, len(tolerance_files_old)):
    if len(tolerance_files_old[iTolFile].split(prefix)) > 2:
        tolerance_files.append(tolerance_files_old[iTolFile].split('_')[1] + '_' + tolerance_files_old[iTolFile].split('_')[2])
    else:
        tolerance_files.append(tolerance_files_old[iTolFile].split(prefix)[1])                                                #

######## 1.PART: create Box Plot
all_averaged_files = []
tolerance_files.insert(6,'')
tolerance_files.insert(13,'')
tolerance_files.insert(20,'')
tolerance_files.insert(27,'')
tolerance_files.insert(34,'')
for iTolerance in range(0, len(tolerance_files)):

    # get empty data in there
    if iTolerance == 6 or iTolerance == 13 or iTolerance == 20 or iTolerance == 27 or iTolerance == 34:
        total_data.append([])

    else:
        # open next .tsv file
        next_tsv = pd.read_csv(tolerance_path + '/' + prefix + tolerance_files[iTolerance], sep='\t')                          #

        # get new .tsv file
        next_tsv = averaging(next_tsv)
        all_averaged_files.append(next_tsv)

        normed_list = []
        for iFile in range(0, len(main_tsv['id'])):
            main_intern = main_tsv['t_intern_ms'][iFile]
            next_intern = next_tsv['t_intern_ms'][iFile]

            # norm all internal + external time by 06_06
            if main_intern == 0:
                quotient = 0
            else:
                quotient = next_intern/main_intern

            # leave out value iff zero
            if quotient == 0:
                'No 0 values allowed!'
            else:
                normed_list.append(quotient)

        # add list to total_data
        total_data.append(normed_list)


# create box_plot
linestyle = (0,(2,5,2,5))
linewidth = 0.1

fontsize = 22 - 4
labelsize = 12
titlesize = 30 + 4

rotation = 45
ax1.set_yscale('log')
print(np.shape(total_data))
bp = ax1.boxplot(total_data, sym='+', widths=0.5, patch_artist=True, positions=range(1,42))#, manage_ticks=True)

####### set more options
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
#ax2.spines['top'].set_linestyle(linestyle)
#ax2.spines['top'].set_linewidth(linewidth)
#ax2.spines['right'].set_linestyle(linestyle)
#ax2.spines['right'].set_linewidth(linewidth)

ax1.set_ylim([0.1,100])
# change colour for each set
color1 = '#66c2a5'
color2 = '#fc8d62'
color3 = '#8da0cb'
color4 = '#e78ac3'
color5 = '#a6d854'
color6 = '#ffd92f'

colors = [color1, color1, color1, color1, color1, color1,
          'white',
          color2, color2, color2, color2, color2, color2,
          'white',
          color3, color3, color3, color3, color3, color3,
          'white',
          color4, color4, color4, color4, color4, color4,
          'white',
          color5, color5, color5, color5, color5, color5,
          'white',
          color6, color6, color6, color6, color6, color6]

#for bplot in bp:
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=2)
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=2)
for median in bp['medians']:
    median.set(color='black', linewidth=2)
for flier in bp['fliers']:
    flier.set(marker='+', color='#e7298a', alpha=0.5)


#ax1.set_title('Comparison of percentiles and median', fontsize=titlesize, fontweight='bold')
ax1.set_ylabel('Relative simulation time', fontsize=fontsize)
ax1.tick_params(labelsize=labelsize)
ax1.set_xticklabels([])
ax1.set_xlim([0,42])
specific_xticks = ax1.xaxis.get_major_ticks()
for iTick in [0,7,14,21,28,35]:
    specific_xticks[iTick].set_visible(False)

# add grit
#ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)


######### 2.PART: add bar plot with failure rate
fontsize = 22 - 4
labelsize = 18 - 5
titlesize = 30

all_percentages = []
counter = 0
for iTolerance in range(0, len(all_averaged_files) + 5):
    if iTolerance in [6, 13, 20, 27, 34]:
        all_percentages.append(0)
        counter += 1
        continue

    # get new .tsv file
    next_tsv = all_averaged_files[iTolerance - counter]
    #next_tsv = averaging(next_tsv)

    # store non-zero and zero values
    non_zero_value_counter = 0
    zero_value_counter = 0
    for iFile in range(0, len(next_tsv['id'])):
        next_intern = next_tsv['t_intern_ms'][iFile]
        if next_intern == 0:
            zero_value_counter += 1
        else:
            non_zero_value_counter += 1

    # store percentage
    all_percentages.append(round(non_zero_value_counter / (non_zero_value_counter + zero_value_counter) * 100, 4))
print(np.shape(all_percentages))

# create bar plot
# colors to use:  '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854'
plot_barplot = ax2.bar(x=list(range(0,41)), height=all_percentages, width=0.5, edgecolor='black', color=colors)

# make top and right boxlines invisible
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
#ax2.spines['top'].set_linestyle(linestyle)
#ax2.spines['top'].set_linewidth(linewidth)
#ax2.spines['right'].set_linestyle(linestyle)
#ax2.spines['right'].set_linewidth(linewidth)

#ax2.set_xscale('log')
ax2.set_xlim([-1,41])
ax2.set_ylim([0,100])
ax2.tick_params(labelsize=labelsize)
ax2.set_ylabel('Success rates [%]', fontsize=labelsize + 5)

# create major and minor ticklabels
#ax2.minorticks_on()
Abs_xTickLabels = [r'$10^{-6}$', '', '', '', '', '', '',
                   r'$10^{-8}$', '', '', '', '', '', '',
                   r'$10^{-10}$', '', '', '', '', '', '',
                   r'$10^{-12}$', '', '', '', '', '', '',
                   r'$10^{-14}$', '', '', '', '', '', '',
                   r'$10^{-16}$']
Rel_xTckLabels = [r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$', '',
                  r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$', '',
                  r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$', '',
                  r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$', '',
                  r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$', '',
                  r'$10^{-6}$', r'$10^{-8}$', r'$10^{-10}$', r'$10^{-12}$', r'$10^{-14}$', r'$10^{-16}$']

ax1.set_xticks(list(range(42)))
ax2.set_xticks(list(range(42)))
minor_list_1 = [x + 0.001 for x in list(range(42))]
ax2.set_xticks(minor_list_1, minor=True)
ax2.set_xticklabels(Abs_xTickLabels, fontsize=labelsize, rotation=rotation)
ax2.set_xticklabels(Rel_xTckLabels, minor=True, fontsize=labelsize, rotation=rotation)
ax2.tick_params(axis='x', which='major', pad=10)
ax2.tick_params(axis='x', which='minor', pad=40)
ax2.text(-0.05, -0.12, 'Abs.tol.: ', fontsize=labelsize, transform=ax2.transAxes)
ax2.text(-0.05, -0.25, 'Rel.tol.: ', fontsize=labelsize, transform=ax2.transAxes)
#ax2.text(-0.05, -0.1, 'Lin. sol.: ', fontsize=fontsize, transform=ax2.transAxes)
#ax2.text(-0.05, -0.18, 'Abs. tol.: ', fontsize=fontsize, transform=ax2.transAxes)
#ax2.text(-0.05, -0.27, 'Rel. tol.: ', fontsize=fontsize, transform=ax2.transAxes)
specific_xticks_major = ax2.xaxis.get_major_ticks()
for iTick in [6,13,20,27,34,41]:
    specific_xticks_major[iTick].set_visible(False)
specific_xticks_minor = ax2.xaxis.get_minor_ticks()
for iTick in [6,13,20,27,34,41]:
    specific_xticks_minor[iTick].set_visible(False)


## better layout
plt.tight_layout()

# change plotting size
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

# save figure
#plt.savefig('../paper_SolverSettings/Figures/Study_2/Tolerances_BoxPlot_BarPlot_' + Multistep_Method + '.pdf')

# show figure
plt.show()