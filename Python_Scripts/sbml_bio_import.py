import os
import urllib.request
import shutil
import sys
import logging
from logging_util import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)

BASE_URL = "https://www.ebi.ac.uk"
BASE_FOLDER = "../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models"


def download_specific_sbml_biomodels_from_jws(model_ids, model_names, base_folder=BASE_FOLDER):
    """
    Download all sbml models to xml files.
    """
    # download every single sbml model
    for model_id in range(0, len(model_ids)):
        # create an own folder for the model
        if not os.path.exists(base_folder + "/" + model_ids[model_id]):
            os.makedirs(base_folder + "/" + model_names[model_id])
        if not os.path.exists(base_folder + "/" + model_names[model_id] + "/sbml_models"):
            os.makedirs(base_folder + "/" + model_names[model_id] + "/sbml_models")
        sbml_file = base_folder + "/" + model_names[model_id] + "/sbml_models/" + model_names[model_id] + ".xml"
        # url to download sbml
        sbml_url = BASE_URL + "/biomodels-main/download?mid=" + model_ids[model_id]
        # download sbml model
        download_sbml_model(sbml_url, sbml_file)


def download_sbml_model(sbml_url, sbml_file):
    """
    Download one sbml model from `sbml_url` to `sbml_file`.
    """
    logger.info(f"  Downloading sbml model from {sbml_url} to file "
                f"{sbml_file}.")

    try:
        with urllib.request.urlopen(sbml_url) as response, \
                open(sbml_file, 'wb') as f:
            shutil.copyfileobj(response, f)
    except Exception as e:
        logger.warn(f"Failed to download sbml model from {sbml_url}, {e}.")


def add_Froehlich2018():
    #basic_path_to_Froehlich2018 = '../Models/all_models'
    #all_models = sorted(os.listdir(basic_path_to_Froehlich2018))
    #if not 'Froehlich2018' in all_models:
    #    print('The name of the repository might have been altered or the folder is missing!')
    #    sys.exit()
    #if not os.path.isfile(basic_path_to_Froehlich2018 + '/Froehlich2018/sbml_models/Froehlich2018.xml'):
    #    print('The name of the .xml file might have been altered or it is missing!')
    #    sys.exit()
    if not os.path.exists(BASE_FOLDER + '/Froehlich2018/sbml_models'):
        os.makedirs(BASE_FOLDER + '/Froehlich2018/sbml_models')
    froehlich2018_url = 'https://raw.githubusercontent.com/ICB-DCM/CS_Signalling_ERBB_RAS_AKT/master/FroehlichKes2018/PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml'
    try:
        with urllib.request.urlopen(froehlich2018_url) as response, \
                open(BASE_FOLDER + '/Froehlich2018/sbml_models/Froehlich2018.xml', 'wb') as f:
            shutil.copyfileobj(response, f)
    except Exception as e:
        logger.warn(f"Failed to download sbml model from {froehlich2018_url}, {e}.")


def add_BioModels_Folder(model_names):
    if not os.path.exists('../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/BioModelsDatabase_models'):
        os.makedirs('../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/BioModelsDatabase_models')
    all_models = sorted(os.listdir(BASE_FOLDER))
    for model in all_models:
        if model in model_names or model == 'Froehlich2018':
            shutil.copytree(BASE_FOLDER + '/' + model, '../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/BioModelsDatabase_models/' + model)