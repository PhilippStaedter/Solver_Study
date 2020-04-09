import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np


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
    adams_models = []
    bdf_models = []
    counters_Adams = []
    counters_BDF = []
    for Multistep in ['Adams', 'BDF']:
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
                        if Multistep == 'Adams':
                            adams_models.append(sedml_model + '_' + sbml_model)
                        elif Multistep == 'BDF':
                            bdf_models.append(sedml_model + '_' + sbml_model)
                    total_counter += 1
            print(tol_exp, counter, total_counter)
            tol_exps.append(tol_exp)
            if Multistep == 'Adams':
                counters_Adams.append(counter)
            elif Multistep == 'BDF':
                counters_BDF.append(counter)
    if iTolerance == 0:
        counter_Tol_0.append(counters_Adams)
        counter_Tol_0.append(counters_BDF)
    elif iTolerance == 1:
        counter_Tol_1.append(counters_Adams)
        counter_Tol_1.append(counters_BDF)
    elif iTolerance == 2:
        counter_Tol_2.append(counters_Adams)
        counter_Tol_2.append(counters_BDF)
    elif iTolerance == 3:
        counter_Tol_3.append(counters_Adams)
        counter_Tol_3.append(counters_BDF)
    elif iTolerance == 4:
        counter_Tol_4.append(counters_Adams)
        counter_Tol_4.append(counters_BDF)


ax = plt.axes()
ax.plot(tol_exps, counter_Tol_0[0], '-*', c='#fdae61', label=f'AM & {all_abs_tol[0]} & {all_rel_tol[0]}')
ax.plot(tol_exps, counter_Tol_0[1], '-*', c='#abd9e9', label=f'BDF & {all_abs_tol[0]} & {all_rel_tol[0]}')
ax.plot(tol_exps, counter_Tol_1[0], '-*', c='#f46d43', label=f'AM & {all_abs_tol[1]} & {all_rel_tol[1]}')
ax.plot(tol_exps, counter_Tol_1[1], '-*', c='#74add1', label=f'BDF & {all_abs_tol[1]} & {all_rel_tol[1]}')
ax.plot(tol_exps, counter_Tol_2[0], '-*', c='#d73027', label=f'AM & {all_abs_tol[2]} & {all_rel_tol[2]}')
ax.plot(tol_exps, counter_Tol_2[1], '-*', c='#4575b4', label=f'BDF & {all_abs_tol[2]} & {all_rel_tol[2]}')
ax.plot(tol_exps, counter_Tol_3[0], '-*', c='#a50026', label=f'AM & {all_abs_tol[3]} & {all_rel_tol[3]}')
ax.plot(tol_exps, counter_Tol_3[1], '-*', c='#313695', label=f'BDF & {all_abs_tol[3]} & {all_rel_tol[3]}')
ax.plot(tol_exps, counter_Tol_4[0], '-*', c='#800000', label=f'AM & {all_abs_tol[4]} & {all_rel_tol[4]}')
ax.plot(tol_exps, counter_Tol_4[1], '-*', c='#2E2D66', label=f'BDF & {all_abs_tol[4]} & {all_rel_tol[4]}')

# local properties
fontsize = 13
ax.set_ylim([-20,400])
ax.set_xticklabels(np.array(['', r'$10^{-20}$', r'$10^{-15}$', r'$10^{-10}$', r'$10^{-5}$', r'$10^{0}$', r'$10^{5}$', r'$10^{10}$']))

# make top and right boxlines invisible
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


plt.legend(loc=2)
#plt.gca().set_title('Comparison of all State Trajectories to JWS for AM vs BDF', fontsize=fontsize)
plt.gca().set_xlabel('Acceptance Threshold for matching State Trajectories', fontsize=fontsize)
plt.gca().set_ylabel('Matching models', fontsize=fontsize)
plt.gcf().tight_layout()

# save plot
#plt.savefig(f"figures_paper/Study_1/compareStateTrajectories_summary_{abs_tol}_{rel_tol}.pdf")

# show plot
plt.show()