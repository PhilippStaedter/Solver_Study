# scatter plot of the prediction settings against amici default settings

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from averageTime import *
from matplotlib.colors import LinearSegmentedColormap

# important paths
base_path = '../paper_SolverSettings/WholeStudy'
Adams_base_path = base_path
BDF_base_path = base_path
#Adams_base_path = '../bachelor_thesis/SolverAlgorithm/Adams'
#BDF_base_path = '../bachelor_thesis/SolverAlgorithm/BDF'

# open .tsv files
list_directory_general = sorted(os.listdir(base_path))
list_directory_adams = []
list_directory_bdf = []
for iFile in range(0, int(len(list_directory_general)/2)):
    list_directory_adams.append(list_directory_general[iFile])
    list_directory_bdf.append(list_directory_general[iFile + int(len(list_directory_general)/2)])
#list_directory_adams = sorted(os.listdir(Adams_base_path))
#list_directory_bdf = sorted(os.listdir(BDF_base_path))

# create list of doubles for scatter plot
adams_bdf_x = []          # red
adams_bdf_y = []
bdf_adams_x = []          # green
bdf_adams_y = []
equal_x = []                    # yellow
equal_y = []
adams_zero_x = []
adams_zero_y = []
bdf_zero_x = []
bdf_zero_y = []
equal_zero_x = []
equal_zero_y = []

for iTsvFile in range(0, len(list_directory_adams)):
    adams_tsv_file = pd.read_csv(Adams_base_path + '/' + list_directory_adams[iTsvFile], sep='\t')
    bdf_tsv_file = pd.read_csv(BDF_base_path + '/' + list_directory_bdf[iTsvFile], sep='\t')

    # average from 210 to 166 models
    adams_tsv_file = averaging(adams_tsv_file)
    bdf_tsv_file = averaging(bdf_tsv_file)

    for iModel in range(0, len(adams_tsv_file['t_intern_ms'])):
        x_adams_data = adams_tsv_file['t_intern_ms'][iModel]
        y_bdf_data = bdf_tsv_file['t_intern_ms'][iModel]

        if x_adams_data != 0 and y_bdf_data != 0:
            if x_adams_data > y_bdf_data:
                bdf_adams_x.append(x_adams_data)
                bdf_adams_y.append(y_bdf_data)
            elif y_bdf_data > x_adams_data:
                adams_bdf_x.append(x_adams_data)
                adams_bdf_y.append(y_bdf_data)
            elif x_adams_data == y_bdf_data:
                equal_x.append(x_adams_data)
                equal_y.append(y_bdf_data)
        elif x_adams_data == 0 and y_bdf_data != 0:
            adams_zero_x.append(45000)
            adams_zero_y.append(y_bdf_data)
        elif x_adams_data != 0 and y_bdf_data == 0:
            bdf_zero_x.append(x_adams_data)
            bdf_zero_y.append(45000)
        elif x_adams_data == 0 and y_bdf_data == 0:
            equal_zero_x.append(45000)
            equal_zero_y.append(45000)

    # display progress
    print(iTsvFile)

# look for the biggest/smallest values
print('adams_bdf_x: ' + str(sorted(adams_bdf_x)[0]))                                            #, reverse=True
print('adams_bdf_y: ' + str(sorted(adams_bdf_y)[0]))
print('bdf_adams_x: ' + str(sorted(bdf_adams_x)[0]))
print('bdf_adams_y: ' + str(sorted(bdf_adams_y)[0]))
#print('equal_x: ' + str(sorted(equal_x, reverse=True)[0]))
#print('equal_y: ' + str(sorted(equal_y, reverse=True)[0]))
print('adams_zero_x: ' + str(sorted(adams_zero_x)[0]))
print('adams_zero_y: ' + str(sorted(adams_zero_y)[0]))
print('bdf_zero_x: ' + str(sorted(bdf_zero_x)[0]))
print('bdf_zero_y: ' + str(sorted(bdf_zero_y)[0]))
print('equal_zero_x: ' + str(sorted(equal_zero_x)[0]))
print('equal_zero_y: ' + str(sorted(equal_zero_y)[0]))
print('len(equal_zero_x): ' + str(len(equal_zero_x)))
print('len(equal_zero_y): ' + str(len(equal_zero_y)))

# plot a scatter plot + diagonal line
linestyle = (0, (2, 5, 2, 5))
linewidth = 1

fontsize = 17
labelsize = 10
titlesize = 30

alpha = 1
marker_size = 2


# create custom colormap
colors_1 = [(1, 0.9, 0.6, 0.3), (0.96, 0.41, 0, 1)]                 # from (red, green, blue, alpha) to (red, green, blue, alpha)
cm_1 = LinearSegmentedColormap.from_list('test', colors_1, N=30)    # number of different color steps
colors_2 = [(0.8, 0.95, 1, 0.3), (0, 0.21, 0.46, 1)]
cm_2 = LinearSegmentedColormap.from_list('test', colors_2, N=30)

# Calculate the point density
# orange
grid_orange = np.vstack([np.log(adams_bdf_x), np.log(adams_bdf_y)])
kde_orange = gaussian_kde(grid_orange)(grid_orange)
# blue
grid_blue = np.vstack([np.log(bdf_adams_x), np.log(bdf_adams_y)])
kde_blue = gaussian_kde(grid_blue)(grid_blue)
# orange edge
grid_orange_edge = np.vstack([np.log(adams_zero_x), np.log(adams_zero_y)])
kde_orange_edge = gaussian_kde(grid_orange_edge)(grid_orange_edge)
# blue edge
grid_blue_edge = np.vstack([np.log(bdf_zero_x), np.log(bdf_zero_y)])
kde_blue_edge = gaussian_kde(grid_blue_edge)(grid_blue_edge)

# Sort the points by density, so that the densest points are plotted last
# orange
ids_orange = kde_orange.argsort()
adams_bdf_x, adams_bdf_y, kde_orange = np.array(adams_bdf_x)[ids_orange], np.array(adams_bdf_y)[ids_orange], np.array(kde_orange)[ids_orange]
#or: adams_bdf_x, adams_bdf_y, kde_orange = [adams_bdf_x[i] for i in ids_orange], [adams_bdf_y[i] for i in ids_orange], [kde_orange for i in ids_orange]
# blue
ids_blue = kde_blue.argsort()
bdf_adams_x, bdf_adams_y, kde_blue = np.array(bdf_adams_x)[ids_blue], np.array(bdf_adams_y)[ids_blue], np.array(kde_blue)[ids_blue]
# orange
ids_orange_edge = kde_orange_edge.argsort()
adams_zero_x, adams_zero_y, kde_orange_edge = np.array(adams_zero_x)[ids_orange_edge], np.array(adams_zero_y)[ids_orange_edge], np.array(kde_orange_edge)[ids_orange_edge]
# blue
ids_blue_edge = kde_blue_edge.argsort()
bdf_zero_x, bdf_zero_y, kde_blue_edge = np.array(bdf_zero_x)[ids_blue_edge], np.array(bdf_zero_y)[ids_blue_edge], np.array(kde_blue_edge)[ids_blue_edge]


# plot scatter plot
ax = plt.axes([0.1, 0.12, 0.8, 0.8])
plt.gcf().subplots_adjust(bottom=0.2)
z = range(0,45000)
plt1 = ax.scatter(adams_bdf_x, adams_bdf_y, s=marker_size, c=kde_orange, cmap=cm_1, label='AM faster: ' + str(round(len(adams_bdf_x) / len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', zorder=10, clip_on=False, alpha=alpha)
plt2 = ax.scatter(bdf_adams_x, bdf_adams_y, s=marker_size, c=kde_blue, cmap=cm_2, label='BDF faster: ' + str(round(len(bdf_adams_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', zorder=10, clip_on=False, alpha=alpha)
plt3 = ax.scatter(equal_x, equal_y, s=marker_size, c='grey', zorder=10, clip_on=False, alpha=alpha) # label='Both are equally good: ' + str(round(len(equal_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %',
plt4 = ax.scatter(adams_zero_x, adams_zero_y, c=kde_orange_edge, cmap=cm_2, marker='D', s=marker_size, facecolors='none', edgecolors='blue', zorder=10, clip_on=False) # label='Adams-Moulton failed to integrate the model: ' + str(round(len(adams_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %',
plt5 = ax.scatter(bdf_zero_x, bdf_zero_y, c=kde_blue_edge, cmap=cm_1, s=marker_size, facecolors='none', edgecolors='orange', marker='D', zorder=10, clip_on=False) # label='BDF failed to integrate the model: ' + str(round(len(bdf_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %',
plt6 = ax.scatter(equal_zero_x, equal_zero_y, s=marker_size, facecolors='none', edgecolors='grey', marker='D', zorder=10, clip_on=False) # label='Both failed to integrate the model: ' + str(round(len(equal_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %',
ax.plot(z, c='black', zorder=20)
ax.set_xlim([0.2, 45000])
ax.set_ylim([0.2, 45000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('AM simulation time [ms]', fontsize=fontsize)
ax.set_ylabel('BDF simulation time [ms]', fontsize=fontsize)
#ax.set_title('Adams-Moulton vs. BDF settings', fontsize=titlesize, fontweight='bold', pad=40)
#lg = ax.legend(loc=2, fontsize=labelsize)
#fr = lg.get_frame()
#fr.set_lw(0.2)

# plot legend manually
ax.plot(0.4, 25000, 'o', fillstyle='full', c='orange', markersize=marker_size)
ax.plot(0.4, 12500, 'o', fillstyle='full', c='blue', markersize=marker_size)
plt.text(0.6, 20000, 'AM faster: ' + str(round(len(adams_bdf_x) / len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', fontsize=labelsize)
plt.text(0.6, 10000, 'BDF faster: ' + str(round(len(bdf_adams_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', fontsize=labelsize)

plt.tick_params(labelsize=labelsize)
plt.gca().set_aspect('equal', adjustable='box')
ax.spines['top'].set_linestyle(linestyle)
ax.spines['top'].set_linewidth(linewidth)
ax.spines['right'].set_linestyle(linestyle)
ax.spines['right'].set_linewidth(linewidth)
ax.spines['top'].set_color('red')
ax.spines['right'].set_color('red')

# write text over axis
ax.text(55000, 1,'only AM failed: ' + str(round(len(adams_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', fontsize=fontsize, rotation=-90)
ax.text(1, 55000, 'only BDF failed: ' + str(round(len(bdf_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', fontsize=fontsize)
ax.text(70000, 4000, 'Both failed:', fontsize=fontsize, rotation=-45)
ax.text(50000, 6000, str(round(len(equal_zero_x)/len(adams_tsv_file['t_intern_ms'])*10/7, 2)) + ' %', fontsize=fontsize, rotation=-45)
#ax.text(0.1, 70000, 'Adams-Moulton vs. BDF settings', fontsize=titlesize, fontweight='bold', pad=30)


# change plotting size
plt.gcf().subplots_adjust(bottom=0.2)
#fig = plt.gcf()
#fig.set_size_inches(18.5, 10.5)
#plt.gcf().subplots_adjust(bottom=0.5)

# better layout
#plt.tight_layout()

# save figure
#plt.savefig('../bachelor_thesis/New_Figures/Figures_study_5/Adams_vs_BDF_2_166SBML.pdf')

# show figure
plt.show()