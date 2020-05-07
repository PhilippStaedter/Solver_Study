# load every working model from the model collection into workspace

import sys
import os
import importlib


def load_specific_model(model_name, explicit_model, skip_indicator):

    # path to one specific model
    if skip_indicator == 0:
        path = '../Models/amici_import'
    elif skip_indicator == 1:
        path = '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sbml2amici/correct_amici_models_paper'
        model_output_dir = path + '/' + model_name + '/' + explicit_model

    # load specific model
    sys.path.insert(0, os.path.abspath(model_output_dir))
    model_module = importlib.import_module(explicit_model)
    model = model_module.getModel()

    return model
