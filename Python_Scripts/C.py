"""Constant definitions."""

###############################################################################
# Base folder

# The output of the study has been conserved in `DIR_BASE`.
DIR_BASE = '../'
# The simulation scripts write their output to `DIR_OUT_BASE` instead.
DIR_OUT_BASE = '../../Benchmarking_of_numerical_ODE_integration_methods/'
# To directly use that output, (un-)comment the below line.
# DIR_BASE = DIR_OUT_BASE


###############################################################################
# Data directories


def dir_data(base=DIR_BASE):
    return base + 'Data/'


def dir_data_wholestudy(base=DIR_BASE):
    return dir_data(base) + 'WholeStudy/'


def dir_data_tolerances(base=DIR_BASE):
    return dir_data(base) + 'TolerancesStudy/'


###############################################################################
# Model directory


def dir_models(base=DIR_BASE):
    return base + 'Models/'


def dir_models_all(base=DIR_BASE):
    return dir_models(base) + 'all_models/'


###############################################################################
# JSON files


def dir_jsonfiles(base=DIR_BASE):
    return base + 'json_files/'
