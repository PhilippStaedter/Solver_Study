# run all models with defined settings

from execute_loadModels import *
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import libsbml
import time
import statistics
import pandas as pd


def simulate(atol, rtol, linSol, iter, solAlg, maxStep, study_indicator, skip_indicator):

    # Create new folder structure for study
    if study_indicator == 1:
        tolerance_path = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/Data/WholeStudy'
    elif study_indicator == 2:
        tolerance_path = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/Data/TolerancesStudy'
    if not os.path.exists(tolerance_path):
        os.makedirs(tolerance_path)

    # save name
    s_atol = str(atol).split('-')[1]
    s_rtol = str(rtol).split('-')[1]
    s_linSol = str(linSol)
    s_iter = str(iter)
    s_solAlg = str(solAlg)

    # create .tsv file
    tsv_table = pd.DataFrame(columns=['id', 't_intern_ms', 't_extern_ms', 'state_variables', 'reactions', 'parameters', 'status', 'error_message'])

    # set row counter for .tsv file
    counter = 0

    # set number of repetitions for simulation
    sim_rep = 40

    # insert specific model properties as strings, e.g.:
    if skip_indicator == 0:
        base_path_sbml2amici = '../Models/amici_import'
        base_path_sedml = '../Models/all_models'
    elif skip_indicator == 0.33:
        base_path_sbml2amici = '../../Models/amici_import'
        base_path_sedml = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models'
    elif skip_indicator == 0.67:
        base_path_sbml2amici = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sbml2amici/correct_amici_models_paper'
        base_path_sedml = '../Models/all_models'
    elif skip_indicator == 1:
        base_path_sbml2amici = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sbml2amici/correct_amici_models_paper'
        base_path_sedml = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models'

    # list of all directories + SBML files
    list_directory_sedml = os.listdir(base_path_sbml2amici)   #(base_path_sedml)
    list_directory_sedml = sorted(list_directory_sedml)
    #list_directory_sedml = list_directory_sedml[27:]
    #list_directory_sedml.remove('ReadMe.MD')

    # list only specific models ---- should only simulate those models where sbml to amicic worked!
    for iModel in list_directory_sedml:

        #iModel = 'Becker_Science2010'

        if os.path.exists(base_path_sbml2amici + '/' + iModel):
            list_files = os.listdir(base_path_sbml2amici + '/' + iModel)
            list_files = sorted(list_files)                                                      # sorted() could maybe change the order needed for later

            #list_directory_xml = [filename for filename in sorted(os.listdir(benchmark_path + '/' + iModel)) if filename.startswith('model_')]
            for iFile in list_files:                                                     #list_directory_xml:

                # Append additional row in .tsv file
                tsv_table = tsv_table.append({}, ignore_index=True)
                #error_file = error_file.append({}, ignore_index=True)

                # save id in .tsv
                tsv_table.loc[counter].id = '{' + iModel + '}' + '_' + '{' + iFile + '}'
                tsv_table.loc[counter].setting = s_solAlg + '_' + s_atol + '_' + s_rtol
                #error_file.loc[counter].id = '{' + iModel + '}' + '_' + '{' + iFile + '}'

                try:
                    # read in SBML file for reactions since AMICI has no functions to count all reactions
                    if os.path.isfile(base_path_sedml + '/' + iModel + '/sbml_models/' + iFile + '.sbml'):
                        file = libsbml.readSBML(base_path_sedml + '/' + iModel + '/sbml_models/' + iFile + '.sbml')
                    else:
                        file = libsbml.readSBML(base_path_sedml + '/' + iModel + '/sbml_models/' + iFile + '.xml')
                    all_properties = file.getModel()
                    num_reactions = all_properties.getNumReactions()
                    tsv_table.loc[counter].reactions = num_reactions

                    # read in model
                    model = all_settings(iModel, iFile, skip_indicator)                   ##################### call function from 'execute_loadModels.py'

                    # save state_variables, reactions and parameters
                    num_states = len(model.getStateNames())
                    num_reactions = len(model.getListOfReactions())
                    num_par = len(model.getParameters())
                    tsv_table.loc[counter].state_variables = num_states
                    tsv_table.loc[counter].parameters = num_par

                    # Create solver instance
                    solver = model.getSolver()

                    # set all settings
                    solver.setAbsoluteTolerance(atol)
                    solver.setRelativeTolerance(rtol)
                    solver.setLinearSolver(linSol)
                    solver.setNonlinearSolverIteration(iter)
                    solver.setLinearMultistepMethod(solAlg)
                    solver.setMaxSteps(maxStep)

                    # clock simulation time while running the simulation using pre-defined settings
                    built_in_time = []
                    external_time = []
                    ind_time = []
                    end_time = []

                    try:
                        # set stability limit detection
                        solver.setStabilityLimitFlag(False)
                        for iSim in range(0, sim_rep):
                            start_time = time.time()

                            sim_data = amici.runAmiciSimulation(model, solver)

                            end_time.append(time.time())                             # x1000 for milliseconds
                            ind_time.append(sim_data['cpu_time'])

                            external_time.append(1000*(end_time[iSim] - start_time))
                            if iSim == 0:
                                built_in_time.append(ind_time[iSim])                                # internal data
                            else:
                                built_in_time.append(ind_time[iSim] - ind_time[iSim - 1])

                        # take status + median of time_vector
                        sim_status = sim_data['status']
                        internal = statistics.median(built_in_time)                                   # median internal data
                        external = statistics.median(external_time)

                        # save status + time data in .tsv
                        tsv_table.loc[counter].status = sim_status
                        tsv_table.loc[counter].t_intern_ms = internal                                    # add internal to .tsv file
                        tsv_table.loc[counter].t_extern_ms = external
                        #error_file.loc[counter].error = sim_status

                        # np.set_printoptions(threshold=8, edgeitems=2)
                        #for key, value in sim_data.items():
                        #    print('%12s: ' % key, value)

                        # plot sim_data
                        # amici.plotting.plotStateTrajectories(sim_data)
                        # amici.plotting.plotObservableTrajectories(sim_data)

                        # save plot in therefore created folder
                        # if not os.path.exists('../sbml2amici/Figures/' + iModel + '/' + iModel):
                        #   os.makedirs('../sbml2amici/Figures/' + iModel + '/' + iModel)
                        # plt.savefig('../sbml2amici/Figures/' + iModel + '/' + iModel + '/' + iFile + '.png')

                        # show plot
                        # plt.show()

                        # raise counter
                        counter = counter + 1

                    except Exception as e:
                        error_info_3 = str(e)
                        # print('Model ' + iModel + '_' + iFile + ' could not be simulated with this setting!')
                        tsv_table.loc[counter].t_intern_ms = 'nan'
                        tsv_table.loc[counter].t_extern_ms = 'nan'
                        tsv_table.loc[counter].error_message = 'Error_3: ' + error_info_3

                        # raise counter
                        counter = counter + 1
                except Exception as e:
                    error_info_2 = str(e)
                    # print('Error_2: Loading Model ' + iModel + '_' + iFile + ' did not work!')
                    tsv_table.loc[counter].t_intern_ms = 'nan'
                    tsv_table.loc[counter].t_extern_ms = 'nan'
                    tsv_table.loc[counter].error_message = 'Error_2: ' + error_info_2

                    # raise counter
                    counter = counter + 1

        ''' (Should not happen, since all models have been simulated before!)
        else:
            'Error_1: Model ' + iModel + ' import to amici did not work!'
        '''
    # save data frame as .tsv file
    tsv_table.to_csv(path_or_buf=tolerance_path + '/' + s_solAlg + '_' + s_iter + '_' + s_linSol + '_' + s_atol + '_' + s_rtol + '.tsv', sep='\t', index=False)

    #return error_file