# execute script loadModels.py + simulate model

from setTime_BioModels import *
from loadModels import *
from changeValues import *
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import libsedml
#from setTime_BenchmarkModels import *


def all_settings(iModel, iFile):

    # insert specific model properties as strings, e.g.:
    base_path_sbml2amici = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sbml2amici/amici_models_newest_version_0.10.19'
    base_path_sedml = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models'
    BioModels_path = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/BioModelsDatabase_models'
    benchmark_path = '../benchmark-models/hackathon_contributions_new_data_format'

    # iFile without extension
    try:
        iFile,_ = iFile.split('.',1)
    except:
        'No extension'

    # run function
    model = load_specific_model(iModel, iFile)                                                          ################ call function from 'loadModels.py'

    if os.path.exists(BioModels_path + '/' + iModel):                                     #(benchmark_path + '/' + iModel):
        sim_start_time, sim_end_time, sim_num_time_points, y_bound = timePointsBioModels(iModel)   #time_array = timePointsBenchmarkModels(iModel, iFile) ################### call function from 'setTime_BioModels.py'
    else:
        # change parameter and species according to SEDML file
        model = changeValues(model, iModel, iFile)                                                      ################# call function from 'changeValues.py'

        # open sedml to get tasks + time courses
        sedml_path = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models/' + iModel + '/' + iModel + '.sedml'

        # tasks
        sedml_file = libsedml.readSedML(sedml_path)

        for iTask in range(0, sedml_file.getNumTasks()):
            all_tasks = sedml_file.getTask(iTask)
            tsk_Id = all_tasks.getId()
            task_name = all_tasks.getName()
            task_modRef = all_tasks.getModelReference()
            task_simReference = all_tasks.getSimulationReference()

            # time courses
            try:
                all_simulations = sedml_file.getSimulation(iTask)
                sim_Id = all_simulations.getId()
            except:                                                                                         # need 'except' clause if more models have same time period
                if all_simulations == None:
                    all_simulations = sedml_file.getSimulation(0)
                    sim_Id = all_simulations.getId()
            try:
                while task_simReference != sim_Id:
                    iTask = iTask + 1
                    all_simulations = sedml_file.getSimulation(iTask)
                    sim_Id = all_simulations.getId()
            except:
                iTask = 0
                while task_simReference != sim_Id:
                    all_simulations = sedml_file.getSimulation(iTask)
                    sim_Id = all_simulations.getId()
                    iTask = iTask + 1

            sim_start_time = all_simulations.getOutputStartTime()
            sim_end_time = all_simulations.getOutputEndTime()
            sim_num_time_points = all_simulations.getNumberOfPoints()

            # load script 'changeValues'
            # model = changeValues(iModel, iFile)

    # set timepoints for which we want to simulate the model
    model.setTimepoints(np.linspace(sim_start_time, sim_end_time, sim_num_time_points))                     #model.setTimepoints(time_array)


    return model