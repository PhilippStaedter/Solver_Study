# script 1 to compare state trajectories from BioModels with trajectories of the simulation created by COPASI
# => creates two important .tsv files

# Attention:    boundary conditions are not being simulated by COPASI!
# all trajectories simulated by COPASI were simulated with absolute and relative tolerances 1e-14 & 1e-14,
# except three models where 1e-12 & 1e-12 had to be used.

from execute_loadModels import *
from JWS_changeValues import *
from colourDataFrame import *
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import libsbml
import libsedml
import time
import statistics
import pandas as pd
import os
import urllib.request
import requests
import json
import itertools


def compStaTraj_BioModels():
    # set settings for simulation
    atol = 1e-3    #1e-1 #1e-3  #1e-6  #1e-12  #1e-6  #1e-14  #1e-16                  # 1e-16
    rtol = 1e-3    #1e-1 #1e-3  #1e-6  #1e-12  #1e-3  #1e-14  #1e-16                  # 1e-8
    linSol = 9
    solAlg = 2                                                                                                              # 1

    # split atol an rtol
    _, save_atol = str('{:.1e}'.format(atol)).split('-')
    _, save_rtol = str('{:.1e}'.format(rtol)).split('-')

    # get all models, but not the README.md if existent
    list_directory_bio = sorted(os.listdir('./StateTrajectories_BioModels_COPASI_Data'))
    if 'README.md' in list_directory_bio:
        list_directory_bio.remove('README.md')
    #list_directory_bio.remove('Froehlich2018')

    for iMod in range(0, len(list_directory_bio)):

        iModel = list_directory_bio[iMod]
        list_files = sorted(os.listdir('./sedml_models/' + iModel + '/sbml_models'))

        for iFile in list_files:

            #iModel = 'Froehlich2018'
            #iFile = iModel + '.xml'

            # iFile without .xml extension
            iFile, extension = iFile.split('.', 1)

            # important paths
            tsv_path = './StateTrajectories_BioModels_COPASI_Data/' + iModel
            save_path = './StateTrajectories_BioModels_COPASI_Data/' + iModel
            sbml_path = './sedml_models/' + iModel + '/sbml_models/' + iFile + '.xml'

            # Open SBML file
            sbml_model = libsbml.readSBML(sbml_path)

            # Get whole model
            model = all_settings(iModel, iFile)

            ######### COPASI simulatiom
            # open new .csv file with COPASI simulation trajectories
            try:
                tsv_file = pd.read_csv(tsv_path + '/Original_COPASI_' + iModel + '_14_14.tsv', sep='\t')
            except:
                tsv_file = pd.read_csv(tsv_path + '/Original_COPASI_' + iModel + '_12_12.tsv', sep='\t')

            # columns names of .tsv file + alter .tsv file
            column_names = list(tsv_file.columns)
            column_names.remove('# Time')
            column_names = column_names[:len(column_names) - 1]
            position = []
            del_counter = 0
            for iName in range(0, len(column_names)):
                if 'Values[' in column_names[iName - del_counter]:
                    column_names.remove(column_names[iName - del_counter])
                    position.append(iName)
                    del_counter += 1
            del tsv_file['# Time']
            tsv_file.drop(tsv_file.columns[len(tsv_file.columns) - 1], axis=1, inplace=True)
            if len(position) != 0:
                tsv_file = tsv_file.drop(tsv_file.columns[position], axis=1)

            ########## model simulation
            # Create solver instance
            solver = model.getSolver()

            # set all settings
            solver.setAbsoluteTolerance(atol)
            solver.setRelativeTolerance(rtol)
            solver.setLinearSolver(linSol)
            solver.setLinearMultistepMethod(solAlg)

            # set stability flag for Adams-Moulton
            if solAlg == 1:
                solver.setStabilityLimitFlag(False)

            # Simulate model
            sim_data = amici.runAmiciSimulation(model, solver)

            # np.set_printoptions(threshold=8, edgeitems=2)
            for key, value in sim_data.items():
                print('%12s: ' % key, value)

            # Get state trajectory
            state_trajectory = sim_data['x']

            # Delete all trajectories for boundary conditions
            delete_counter = 0
            all_properties = sbml_model.getModel()
            for iSpec in range(0, all_properties.getNumSpecies()):
                all_species = all_properties.getSpecies(iSpec)
                if all_species.getBoundaryCondition() == True:
                    state_trajectory = state_trajectory.transpose()
                    if delete_counter == 0:
                        state_trajectory = np.delete(state_trajectory, iSpec, 0)
                    else:
                        state_trajectory = np.delete(state_trajectory, iSpec - delete_counter, 0)
                    state_trajectory = state_trajectory.transpose()
                    delete_counter = delete_counter + 1

            # Adjustment for the 'Froehlich2018' model
            if iModel == 'Froehlich2018':
                state_trajectory = np.delete(state_trajectory, np.s_[1177:1178], axis=1)

            # Convert ndarray 'state-trajectory' to data frame + save it and modified COPASI file
            df_state_trajectory = pd.DataFrame(columns=column_names, data=state_trajectory)
            df_state_trajectory.to_csv(f"{save_path}/AMICI_{iModel}_{save_atol}_{save_rtol}.tsv", sep='\t', index=False)
            tsv_file.to_csv(f"{save_path}/COPASI_{iModel}_{save_atol}_{save_rtol}.tsv", sep='\t', index=False)

            # to know where one is
            print('Model ' + iModel + ' successfully completed!')

# call function with no models to delete
compStaTraj_BioModels()