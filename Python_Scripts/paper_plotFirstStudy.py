# combined plot for study one

import matplotlib.pyplot as plt
import pandas as pd
from averageTime import *
import numpy as np



# all axes objects
fontsize = 9
labelsize = 9

rotation = 90
left = 0.07             #0.07
bottom = 0.75
width = 0.4
height = 0.22
row_factor = 0.45
column_factor = 0.11
rotation_factor = 70
alpha = 1

ax1 = plt.axes([0.57, bottom, width, height])
ax2 = plt.axes([0.57, bottom - column_factor - height, width, height])
ax3 = plt.axes([0.57, bottom - 2 * column_factor - 2 * height, width, height])
ax4 = plt.axes([left, 0.09, 0.41, 0.88]) #0.12

############## plot scatter plot ################
all_log_abs_tol = ['03', '06', '06', '16', '12']
all_log_rel_tol = ['03', '03', '06', '08', '12']
all_abs_tol = [r'$10^{-3}$', r'$10^{-6}$', r'$10^{-6}$', r'$10^{-16}$', r'$10^{-12}$']
all_rel_tol = [r'$10^{-3}$', r'$10^{-3}$', r'$10^{-6}$', r'$10^{-8}$', r'$10^{-12}$']

counter_Tol_0 = []
counter_Tol_1 = []
counter_Tol_2 = []
counter_Tol_3 = []
counter_Tol_4 = []
for iTolerance in range(0, len(all_log_abs_tol)):
    abs_tol = all_log_abs_tol[iTolerance]
    rel_tol = all_log_rel_tol[iTolerance]
    bdf_models = []
    counters_BDF = []
    for Multistep in ['BDF']:
        print(Multistep)
        tol_exps = []
        for tol_exp in range(-20, 10):                # range(-20,10)
            tol = 10**tol_exp
            counter = 0
            total_counter = 0
            tol_str = f"{float(tol):.0e}"
            base_dir = f"json_files_all_results_{Multistep}_{abs_tol}_{rel_tol}/json_files_{tol_str}_{tol_str}"
            sedml_models = os.listdir(base_dir)
            for sedml_model in sedml_models:
                sedml_dir = base_dir + "/" + sedml_model
                sbml_models = os.listdir(sedml_dir)
                for sbml_model in sbml_models:
                    sbml_dir = sedml_dir + "/" + sbml_model
                    filename = sbml_dir + "/" + "whole_error.csv"
                    df = pd.read_csv(filename, sep='\t')
                    if df['trajectories_match'][0] == True:
                        counter += 1
                        bdf_models.append(sedml_model + '_' + sbml_model)
                    total_counter += 1
            print(tol_exp, counter, total_counter)
            tol_exps.append(tol_exp)
            counters_BDF.append(counter)
    if iTolerance == 0:
        counter_Tol_0.append(counters_BDF)
    elif iTolerance == 1:
        counter_Tol_1.append(counters_BDF)
    elif iTolerance == 2:
        counter_Tol_2.append(counters_BDF)
    elif iTolerance == 3:
        counter_Tol_3.append(counters_BDF)
    elif iTolerance == 4:
        counter_Tol_4.append(counters_BDF)


ax4.plot(tol_exps, counter_Tol_0[0], '-*', c='#d7191c', label=f'Abs. tol.: {all_abs_tol[0]}, Rel. tol.: {all_rel_tol[0]}')
ax4.plot(tol_exps, counter_Tol_1[0], '-*', c='#fdae61', label=f'Abs. tol.: {all_abs_tol[1]}, Rel. tol.: {all_rel_tol[1]}')
ax4.plot(tol_exps, counter_Tol_2[0], '-*', c='#abd9e9', label=f'Abs. tol.: {all_abs_tol[2]}, Rel. tol.: {all_rel_tol[2]}')
ax4.plot(tol_exps, counter_Tol_3[0], '-*', c='#2c7bb6', label=f'Abs. tol.: {all_abs_tol[3]}, Rel. tol.: {all_rel_tol[3]}')
ax4.plot(tol_exps, counter_Tol_4[0], '-*', c='#2E2D66', label=f'Abs. tol.: {all_abs_tol[4]}, Rel. tol.: {all_rel_tol[4]}')

# plot text 'B'
ax4.text(-0.16, 1, 'A', fontsize=labelsize + 3, transform=ax4.transAxes)

# local properties
ax4.set_ylim([-20,400])
ax4.set_xticks(np.array([-20, -15, -10, -5, 0, 5, 10]), (r'$10^{-20}$', r'$10^{-15}$', r'$10^{-10}$', r'$10^{-5}$', r'$10^{0}$', r'$10^{5}$', r'$10^{10}$'))

# global properties
ax4.legend(loc=2, fontsize=fontsize - 1)
ax4.tick_params(labelsize=fontsize)

# make top and right boxlines invisible
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)

# more properties
plt.gca().set_xlabel("Acceptance Threshold for matching State Trajectories", fontsize=fontsize)
plt.gca().set_ylabel("Matching models", fontsize=fontsize)
plt.gcf().tight_layout()


################ plot bar plot #################
# open two .tsv file
path = '../sbml2amici/NEW_stat_reac_par_paper.tsv'
tsv_file = pd.read_csv(path, sep='\t')
tsv_file = averaging(tsv_file)

# no BioModels yet + delete nans at the end
tsv_file = tsv_file#[27:423]
tsv_file = tsv_file.reset_index()
del tsv_file['index']


# take number of states for those models that worked
data_states_ok = []
for iLine in range(0, len(tsv_file['id'])):
    data_states_ok.append(np.log10(tsv_file['state_variables'][iLine])   )              # no 'ignore_index=True' permitted

# take number of reactions for those models that worked
data_reactions_ok = []
for iLine in range(0, len(tsv_file['id'])):
    data_reactions_ok.append(np.log10(tsv_file['reactions'][iLine]))

# take number of parameters for those models that worked
data_parameters_ok = []
for iLine in range(0, len(tsv_file['id'])):
    data_parameters_ok.append(np.log10(tsv_file['parameters'][iLine]))


# histogram of states
#bins=np.logspace(0,4,200)
bins = 80

#plt.subplot(3,1,1)
# add title
#plt.title('Basic properties of all models', fontsize=20)
plot1 = ax1.hist(x=data_states_ok, range=[-1,4], bins=bins, log=True) # range=[0,250],
#ax1.set_xscale('log')
ax1.set_xlim((-0.1, 4)) #2000 #250 #100                                                                       # Froehlich2018: 1396
ax1.set_ylim((0.5, 150))
ax1.set_xticklabels(['', r'$10^{0}$', '', r'$10^{1}$', '', r'$10^{2}$', '', r'$10^{3}$', '', r'$10^{4}$'])
ax1.set_xlabel('Number of state variables', fontsize=fontsize)
ax1.set_ylabel('Number of models', fontsize=fontsize)
ax1.tick_params(labelsize=labelsize)

# plot text 'A'
ax1.text(-0.22, 1, 'B', fontsize=labelsize + 3, transform=ax1.transAxes)

# histogram of reactions
#plt.subplot(3,1,2)
plot2 = ax2.hist(x=data_reactions_ok, range=[-1,4], bins=bins, log=True) # range=[0,600],
#ax2.set_xscale('log')
ax2.set_xlim((-0.1, 4)) #3000 #600 #100                                                                       # Froehlich2018: 2686
ax2.set_ylim((0.5, 150))
ax2.set_xticklabels(['', r'$10^{0}$', '', r'$10^{1}$', '', r'$10^{2}$', '', r'$10^{3}$', '', r'$10^{4}$'])
ax2.set_xlabel('Number of reactions', fontsize=fontsize)
ax2.set_ylabel('Number of models', fontsize=fontsize)
ax2.tick_params(labelsize=labelsize)

# histogram of parameters
#plt.subplot(3,1,3)
plot3 = ax3.hist(x=data_parameters_ok, range=[-1,4], bins=bins, log=True) # range=[0,350],
#ax3.set_xscale('log')
ax3.set_xlim((-0.1, 4)) #5000 #350 #1000                                                                       # Froehlich2018: 4704
ax3.set_ylim((0.5, 150))
ax3.set_xticklabels(['', r'$10^{0}$', '', r'$10^{1}$', '', r'$10^{2}$', '', r'$10^{3}$', '', r'$10^{4}$'])
ax3.set_xlabel('Number of parameters',  fontsize=fontsize)
ax3.set_ylabel('Number of models', fontsize=fontsize)
ax3.tick_params(labelsize=labelsize)

# make all top and right boxlines inviisible
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# better layout
plt.tight_layout()

# change plotting size
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)


# better layout
plt.tight_layout()

# change plotting size
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

# save figure
#plt.savefig('../paper_SOlverSettings/Figures/Study_1/stat_reac_par_ylog.pdf')

# show figure
plt.show()