# script 1 to compare state trajectories from JWS with trajectories of the AMICI simulation
# => creates two important .tsv files

# Attention:    boundary conditions are not being simulated by JWS!

# Note: libsedml must be imported before libsbml for whatever reason

from execute_loadModels import all_settings
from JWS_changeValues import JWS_changeValues

import amici.plotting
import numpy as np
import libsbml
import libsedml
import pandas as pd
import os
import urllib.request
import requests
import json
import itertools

from C import (
    DIR_MODELS_SEDML, DIR_MODELS_BIOMODELS, DIR_MODELS_AMICI, DIR_MODELS_JSON)


def compStaTraj(delete_counter):
    # set settings for simulation
    for solAlg in [1, 2]:
        linSol = 9
        if solAlg == 1:
            MultistepMethod = 'Adams'
            Tolerance_combination = [
                [1e-3,1e-3], [1e-6,1e-3], [1e-6,1e-6], [1e-12,1e-12],
                [1e-16,1e-8]]
        elif solAlg == 2:
            MultistepMethod = 'BDF'
            Tolerance_combination = [
                [1e-3,1e-3], [1e-6,1e-3], [1e-6,1e-6], [1e-12,1e-12],
                [1e-16,1e-8], [1e-14,1e-14], [1e-16,1e-16]]

        for iTolerance in Tolerance_combination:
            # split atol and rtol for naming purposes
            atol_exp = str(iTolerance[0])
            rtol_exp = str(iTolerance[1])
            if len(atol_exp) != 2:
                atol_exp = '0' + atol_exp
            if len(rtol_exp) != 2:
                rtol_exp = '0' + rtol_exp

            # get name of jws reference
            url = "https://jjj.bio.vu.nl/rest/models/?format=json"
            view_source = requests.get(url)
            json_string = view_source.text
            json_dictionary = json.loads(json_string)

            # get all models
            list_directory_amici = sorted(os.listdir(DIR_MODELS_AMICI))
            if delete_counter != 0:
                del list_directory_amici[0:delete_counter]

            for iMod in range(0, len(list_directory_amici)):
                iModel = list_directory_amici[iMod]
                list_files = sorted(os.listdir(os.path.join(
                    DIR_MODELS_SEDML, iModel, 'sbml_models')))

                for iFile in list_files:
                    print(iModel, iFile, MultistepMethod, linSol, iTolerance)
                    # iFile without .sbml extension
                    iFile, extension = iFile.split('.', 1)

                    # important paths
                    json_save_path = os.path.join(
                        DIR_MODELS_JSON,
                        f'json_files_{MultistepMethod}_{atol_exp}_{rtol_exp}',
                        iModel, iFile)
                    sedml_path = os.path.join(
                        DIR_MODELS_SEDML, iModel, iModel +'.sedml')
                    sbml_path = os.path.join(
                        DIR_MODELS_SEDML, iModel, 'sbml_models',
                        iFile + '.sbml')
                    BioModels_path = DIR_MODELS_BIOMODELS

                    if os.path.exists(os.path.join(BioModels_path, iModel)):
                        print('Model is not part of JWS-database!')
                    else:
                        # Open SBML file
                        sbml_model = libsbml.readSBML(sbml_path)

                        # get right model reference from sbml model
                        parse_name_model = sbml_model.getModel().getId()
                        for iCount in range(0, len(json_dictionary)):
                            parse_name_jws = json_dictionary[iCount]['slug']
                            if parse_name_model == parse_name_jws:
                                model_reference = json_dictionary[iCount]['slug']
                                break
                        # elements in json_dictionary are only lower case --- the sbml model has upper case models
                        try:
                            model_reference
                        except:
                            wrong_model_name = ["".join(x) for _, x in itertools.groupby(parse_name_model,
                                                                                         key=str.isdigit)]
                            if wrong_model_name[0].islower() == False:
                                correct_model_letters = wrong_model_name[0].lower()
                                correct_model_name = correct_model_letters + wrong_model_name[1]
                                for iCount in range(0, len(json_dictionary)):
                                    parse_name_jws = json_dictionary[iCount]['slug']
                                    if correct_model_name == parse_name_jws:
                                        model_reference = json_dictionary[iCount]['slug']
                                        break
                        # check if 'all_settings' works
                        try:
                            # Get whole model
                            model = all_settings(iModel, iFile)

                            # create folder
                            if not os.path.exists(json_save_path):
                                os.makedirs(json_save_path)
                        except:
                            print('Model ' + iModel + ' extension is missing!')
                            continue

                        ######### jws simulation
                        # Get time data with num_time_points == 100
                        t_data = model.getTimepoints()
                        sim_start_time = t_data[0]
                        sim_end_time = t_data[len(t_data) - 1]
                        sim_num_time_points = 101
                        model.setTimepoints(np.linspace(sim_start_time, sim_end_time, sim_num_time_points))

                        # Open sedml file
                        sedml_model = libsedml.readSedML(sedml_path)

                        # import all changes from SEDML
                        list_of_strings = JWS_changeValues(iFile, sedml_model)

                        # Get Url with all changes
                        # <species 1>=<amount>
                        # <parameter 1>=<value>, compartment == parameter (in this case)
                        url = 'https://jjj.bio.vu.nl/rest/models/' + \
                              model_reference + '/time_evolution?time_end=' + \
                              str(sim_end_time) + ';species=all;'

                        for iStr in list_of_strings:
                            url = url + iStr

                        #### Save .json file
                        json_file = os.path.join(
                            json_save_path, iFile + '_JWS_simulation.json')
                        urllib.request.urlretrieve(url, json_file)

                        #### write as .csv file
                        json_2_csv = pd.read_json(json_file)
                        tsv_file = os.path.join(
                            json_save_path, iFile + '_JWS_simulation.csv')
                        json_2_csv.to_csv(tsv_file, sep='\t', index=False)

                        # open new .csv file
                        df = pd.read_csv(tsv_file, sep='\t')

                        # columns names of .tsv file
                        column_names = list(df.columns)
                        column_names.remove('time')
                        del df['time']

                        ########## model simulation
                        # Create solver instance
                        solver = model.getSolver()

                        # set all settings
                        solver.setAbsoluteTolerance(iTolerance[0])
                        solver.setRelativeTolerance(iTolerance[1])
                        solver.setLinearSolver(linSol)
                        solver.setLinearMultistepMethod(solAlg)

                        # set stability flag for Adams-Moulton
                        if solAlg == 1:
                            solver.setStabilityLimitFlag(False)

                        # Simulate model
                        sim_data = amici.runAmiciSimulation(model, solver)

                        # print some values
                        #for key, value in sim_data.items():
                        #       print('%12s: ' % key, value)

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

                        # Convert ndarray 'state-trajectory' to data frame and save it
                        try:
                            df_state_trajectory = pd.DataFrame(columns=column_names, data=state_trajectory)
                        except:
                            print('Try again for model ' + list_directory_amici[iMod] + '_' + iFile)
                            compStaTraj(iMod)
                        df_state_trajectory.to_csv(os.path.join(
                            json_save_path, iFile + '_model_simulation.csv'),
                            sep='\t')

# call function, starting with no models to delete
compStaTraj(0)