# scatter plot of the prediction settings against amici default settings

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from averageTime import *
from matplotlib.colors import LinearSegmentedColormap

# important paths
base_path_LSODA = '../paper_SolverSettings/WholeStudy_LSODA'
base_path = '../paper_SolverSettings/WholeStudy'
LSODA_base_path = base_path_LSODA
AM_base_path = base_path

# general plotting settings
plt.rcParams['figure.figsize'] = [15.0, 7.0]
plt.rcParams['figure.dpi'] = 80
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['font.size'] = 17

# open .tsv files --- AM
list_directory_general = sorted(os.listdir(base_path))
list_directory_am = []
for iFile in range(0, int(len(list_directory_general)/2)):
    am_split = list_directory_general[iFile + int(len(list_directory_general)/2)].split('_')
    if am_split[0] == '1' and am_split[1] == '2'and am_split[2] == '9':
        list_directory_am.append(list_directory_general[iFile + int(len(list_directory_general)/2)])

# open .tsv files --- LSODA
list_directory_COPASI = sorted(os.listdir(base_path_LSODA))
list_directory_lsoda = []
for iFile in range(0, len(list_directory_COPASI)):
    lsoda_split = list_directory_COPASI[iFile].split('_')
    if lsoda_split[0] == '(1,2)' and lsoda_split[1] == '2'and lsoda_split[2] == '1':
        list_directory_lsoda.append(list_directory_COPASI[iFile])

# create list of doubles for scatter plot
lsoda_am_x = []          # red
lsoda_am_y = []
am_lsoda_x = []          # green
am_lsoda_y = []
ratio_lsoda_am = []
ratio_am_lsoda = []
ratio_equal = []
ratio_lsoda_zero = []
ratio_am_zero = []
ratio_equal_zero = []
num_states_lsoda_am = []
num_states_am_lsoda = []
num_states_equal = []
num_states_lsoda_zero = []
num_states_am_zero = []
num_states_equal_zero = []
equal_x = []              # yellow
equal_y = []
lsoda_zero_x = []
lsoda_zero_y = []
am_zero_x = []
am_zero_y = []
equal_zero_x = []
equal_zero_y = []

for iTsvFile in range(0, len(list_directory_lsoda)):
    lsoda_tsv_file = pd.read_csv(LSODA_base_path + '/' + list_directory_lsoda[iTsvFile], sep='\t')
    am_tsv_file = pd.read_csv(AM_base_path + '/' + list_directory_am[iTsvFile], sep='\t')

    # average from 211 to 167 models
    lsoda_tsv_file = averaging(lsoda_tsv_file)
    am_tsv_file = averaging(am_tsv_file)

    for iModel in range(0, len(lsoda_tsv_file['t_intern_ms'])):
        x_lsoda_data = lsoda_tsv_file['t_intern_ms'][iModel] * 1000
        y_am_data = am_tsv_file['t_intern_ms'][iModel]
        states_data = lsoda_tsv_file['state_variables'][iModel]

        if x_lsoda_data != 0 and y_am_data != 0:
            if x_lsoda_data > y_am_data:
                am_lsoda_x.append(x_lsoda_data)
                am_lsoda_y.append(y_am_data)
                ratio_lsoda_am.append(np.log10(x_lsoda_data / y_am_data))
                num_states_lsoda_am.append(states_data)
            elif y_am_data > x_lsoda_data:
                lsoda_am_x.append(x_lsoda_data)
                lsoda_am_y.append(y_am_data)
                ratio_am_lsoda.append(np.log10(x_lsoda_data / y_am_data))
                num_states_am_lsoda.append(states_data)
            elif x_lsoda_data == y_am_data:
                equal_x.append(x_lsoda_data)
                equal_y.append(y_am_data)
                ratio_equal.append(np.log10(x_lsoda_data / y_am_data))
                num_states_equal.append(states_data)

        elif x_lsoda_data == 0 and y_am_data != 0:
            lsoda_zero_x.append(300000)
            lsoda_zero_y.append(y_am_data)
            ratio_lsoda_zero.append(np.log10(1000))
            num_states_lsoda_zero.append(states_data)

        elif x_lsoda_data != 0 and y_am_data == 0:
            am_zero_x.append(x_lsoda_data)
            am_zero_y.append(300000)
            ratio_am_zero.append(-np.log10(1000))
            num_states_am_zero.append(states_data)

        elif x_lsoda_data == 0 and y_am_data == 0:
            equal_zero_x.append(300000)
            equal_zero_y.append(300000)
            ratio_equal_zero.append(np.log10(float('nan')))
            num_states_equal_zero.append(states_data)

    # display progress
    print(iTsvFile)

# print some interesting properties --- look for the biggest/smallest values
print('lsoda_am_x_smallest: ' + str(sorted(lsoda_am_x)[0]))
print('lsoda_am_y_smallest: ' + str(sorted(lsoda_am_y)[0]))
print('am_lsoda_x_smallest: ' + str(sorted(am_lsoda_x)[0]))
print('am_lsoda_y_smallest: ' + str(sorted(am_lsoda_y)[0]))
print('lsoda_zero_x_smallest: ' + str(sorted(lsoda_zero_x)[0]))
print('lsoda_zero_y_smallest: ' + str(sorted(lsoda_zero_y)[0]))
print('equal_zero_x_smallest: ' + str(sorted(equal_zero_x)[0]))
print('equal_zero_y_smallest: ' + str(sorted(equal_zero_y)[0]))
print('lsoda_am_x_largest: ' + str(sorted(lsoda_am_x, reverse=True)[0]))
print('lsoda_am_y_largest: ' + str(sorted(lsoda_am_y, reverse=True)[0]))
print('am_lsoda_x_largest: ' + str(sorted(am_lsoda_x, reverse=True)[0]))
print('am_lsoda_y_largest: ' + str(sorted(am_lsoda_y, reverse=True)[0]))
print('lsoda_zero_x_largest: ' + str(sorted(lsoda_zero_x, reverse=True)[0]))
print('lsoda_zero_y_largest: ' + str(sorted(lsoda_zero_y, reverse=True)[0]))
print('equal_zero_x_largest: ' + str(sorted(equal_zero_x, reverse=True)[0]))
print('equal_zero_y_largest: ' + str(sorted(equal_zero_y, reverse=True)[0]))
print('len(equal_zero_x): ' + str(len(equal_zero_x)))
print('len(equal_zero_y): ' + str(len(equal_zero_y)))
print('ratio_am_lsoda_largest: ' + str(sorted(ratio_am_lsoda, reverse=True)[0]))
print('ratio_am_lsoda_samllest: ' + str(sorted(ratio_am_lsoda)[0]))
print('ratio_lsoda_am_largest: ' + str(sorted(ratio_lsoda_am, reverse=True)[0]))
print('ratio_lsoda_am_smallest: ' + str(sorted(ratio_lsoda_am)[0]))

# plot a scatter plot + diagonal line
linestyle = (0, (2, 5, 2, 5))
linewidth = 1
fontsize = 17
labelsize = 10
titlesize = 30
alpha = 1
marker_size = 2

# create custom colormap
colors_1 = [(0.83, 0.83, 0.83, 0.3), (0.3, 0.2, 0.25, 1)]
cm_1 = LinearSegmentedColormap.from_list('test', colors_1, N=30)
colors_2 = [(0.8, 0.95, 1, 0.3), (0, 0.21, 0.46, 1)]
cm_2 = LinearSegmentedColormap.from_list('test', colors_2, N=30)

# Calculate the point density
# grey
grid_grey = np.vstack([np.log(lsoda_am_x), np.log(lsoda_am_y)])
kde_grey = gaussian_kde(grid_grey)(grid_grey)
# blue
grid_blue = np.vstack([np.log(am_lsoda_x), np.log(am_lsoda_y)])
kde_blue = gaussian_kde(grid_blue)(grid_blue)

# Sort the points by density, so that the densest points are plotted last
# grey
ids_grey = kde_grey.argsort()
lsoda_am_x, lsoda_am_y, kde_grey = np.array(lsoda_am_x)[ids_grey], np.array(lsoda_am_y)[ids_grey], np.array(kde_grey)[ids_grey]
# blue
ids_blue = kde_blue.argsort()
am_lsoda_x, am_lsoda_y, kde_blue = np.array(am_lsoda_x)[ids_blue], np.array(am_lsoda_y)[ids_blue], np.array(kde_blue)[ids_blue]


# plot scatter plot
ax = plt.axes([0.08, 0.11, 0.37, 0.84])
ax2 = plt.axes([0.57, 0.11, 0.40, 0.84])
z = range(0, 300000)
plt1 = ax.scatter(lsoda_am_x, lsoda_am_y, s=marker_size, c=kde_grey, cmap=cm_1, label='LSODA faster: ' + str(round(len(lsoda_am_x) / len(lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %', zorder=10, clip_on=False, alpha=alpha)
plt2 = ax.scatter(am_lsoda_x, am_lsoda_y, s=marker_size, c=kde_blue, cmap=cm_2, label='AM faster: ' + str(round(len(am_lsoda_x)/len(lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %', zorder=10, clip_on=False, alpha=alpha)
plt3 = ax.scatter(equal_x, equal_y, s=marker_size, c='grey', zorder=10, clip_on=False, alpha=alpha)
plt4 = ax.scatter(lsoda_zero_x, lsoda_zero_y, c='#4c3340', cmap=cm_2, marker='D', s=marker_size, facecolors='none', edgecolors='blue', zorder=10, clip_on=False)
plt5 = ax.scatter(am_zero_x, am_zero_y, c='blue', cmap=cm_1, s=marker_size, facecolors='none', edgecolors='#4c3340', marker='D', zorder=10, clip_on=False)
plt6 = ax.scatter(equal_zero_x, equal_zero_y, s=marker_size, facecolors='none', edgecolors='grey', marker='D', zorder=10, clip_on=False)
ax.plot(z, c='black', zorder=20)
ax.set_xlim([0.1, 300000])
ax.set_ylim([0.1, 300000])
ax.set_xscale('log')
ax.set_yscale('log')
'''
ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
               2, 3, 4, 5, 6, 7, 8, 9, 10,
               20, 30, 40, 50, 60, 70, 80, 90, 100,
               200, 300, 400, 500, 600, 700, 800, 900, 1000,
               2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
               20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000,
               200000, 300000])
'''
ax.set_xlabel('LSODA simulation time [ms]', fontsize=fontsize)
ax.set_ylabel('AM simulation time [ms]', labelpad=-5, fontsize=fontsize)

# plot legend manually
ax.plot(0.4, 100000, 'o', fillstyle='full', c='#4c3340', markersize=marker_size)
ax.plot(0.4, 37000, 'o', fillstyle='full', c='blue', markersize=marker_size)
ax.text(0.6, 80000, 'LSODA faster: ' + str(round(len(lsoda_am_x) / len(lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %', fontsize=fontsize)
ax.text(0.6, 30000, 'AM faster: ' + str(round(len(am_lsoda_x)/len(lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %', fontsize=fontsize)

plt.tick_params(labelsize=labelsize)
ax.spines['top'].set_linestyle(linestyle)
ax.spines['top'].set_linewidth(linewidth)
ax.spines['right'].set_linestyle(linestyle)
ax.spines['right'].set_linewidth(linewidth)
ax.spines['top'].set_color('grey')
ax.spines['right'].set_color('grey')

# write text over axis
ax.text(400000, 150,'only LSODA failed: ' + str(round(len(lsoda_zero_x) / len(
    lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %',
        fontsize=fontsize, rotation=-90, va='center')
ax.text(100, 400000, 'only AM failed: ' + str(round(len(am_zero_x) / len(
    lsoda_tsv_file['t_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %',
        fontsize=fontsize, ha='center')
ax.text(500000, 400000, 'Both failed:', fontsize=fontsize,
        ha='center')
ax.text(1500000, 200000, str(round(len(equal_zero_x) / len(lsoda_tsv_file[
    't_intern_ms'])*100/len(list_directory_lsoda), 2)) + ' %',
        fontsize=fontsize, ha='center', va='center')

# plot text 'C'
ax.text(-0.18, 1, 'C', fontsize=fontsize + 5, transform=ax.transAxes)


# choose second axes object
plt.sca(ax2)

# create and adapt spines
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# plot left and right "spines"
plt.plot([-np.log10(1000), -np.log10(1000)], [0.7, 2000], '--', color='grey',
         linestyle=linestyle, linewidth=linewidth)
plt.plot([np.log10(1000), np.log10(1000)], [0.7, 2000], '--', color='grey',
         linestyle=linestyle, linewidth=linewidth)

# plot data
plt.scatter(ratio_lsoda_am, num_states_lsoda_am,
            s=marker_size, c=kde_blue, cmap=cm_2, alpha=alpha,
            zorder=10, clip_on=False)
plt.scatter(ratio_am_lsoda, num_states_am_lsoda,
            s=marker_size, c=kde_grey, cmap=cm_1, alpha=alpha,
            zorder=10, clip_on=False)
plt.scatter(ratio_equal, num_states_equal,
            s=marker_size, c='grey', alpha=alpha,
            zorder=100, clip_on=False)
plt.scatter(ratio_lsoda_zero, num_states_lsoda_zero,
            s=marker_size, c='blue', cmap=cm_2, alpha=alpha, marker='D',
            zorder=10, clip_on=False)
plt.scatter(ratio_am_zero, num_states_am_zero,
            s=marker_size, c='#4c3340', cmap=cm_1, alpha=alpha, marker='D',
            zorder=10, clip_on=False)
plt.scatter(ratio_equal_zero, num_states_equal_zero,
            s=marker_size, c='grey', cmap=cm_1, alpha=alpha, marker='D',
            facecolors='none', edgecolors='grey', zorder=100, clip_on=False)


# plot central spine
plt.plot([0, 0], [0.7, 2000], 'k-', linewidth=linewidth)
plt.plot([-.05, .05], [1000, 1000], 'k-', linewidth=linewidth)
for iMajor in range(3):
    plt.plot([-.05, .05], [10**iMajor, 10**iMajor], 'k-', linewidth=linewidth)
    for iMinor in range(2, 10):
        plt.plot([-.02, .02], [iMinor * 10**iMajor, iMinor * 10**iMajor],
                 'k-',linewidth=linewidth)

# formatting
ax2.set_yscale('log')
ax2.set_ylim((0.7, 1900))
ax2.set_xlim((-np.log10(1000), np.log10(1000)))
ax2.set_yticks([])
plt.minorticks_off()
ax2.text(0, 2000, 'Number of state variables', fontsize=fontsize, ha='center')
ax2.text(-0.1, 1, '$10^0$', fontsize=fontsize, va='center', ha='right')
ax2.text(-0.1, 10, '$10^1$', fontsize=fontsize, va='center', ha='right')
ax2.text(-0.1, 100, '$10^2$', fontsize=fontsize, va='center', ha='right')
ax2.text(-0.1, 1000, '$10^3$', fontsize=fontsize, va='center', ha='right')
ax2.set_xlabel('Computation time ratio LSODA / AM', fontsize=fontsize)
plt.xticks([-np.log10(1000), -2, -1, 0, 1, 2, np.log10(1000)],
           ['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1', '$10^1$', '$10^2$',
            '$10^3$'], fontsize=fontsize)
ax2.text(-np.log10(1500), 20, 'AM failed', fontsize=fontsize, rotation=90, ha='center')
ax2.text(np.log10(1500), 20, 'LSODA failed', fontsize=fontsize, rotation=-90, ha='center')

# plot text 'D'
ax.text(-0.08, 1, 'D', fontsize=fontsize + 5, transform=ax2.transAxes)


# change plotting size
plt.gcf().subplots_adjust(bottom=0.2)

# better layout
plt.tight_layout()

# show figure
plt.show()
