import libsedml
import os
import urllib.request
import shutil
import json
import logging
from logging_util import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)

BASE_URL = "https://www.ebi.ac.uk"
BASE_FOLDER = "../../Assessment_of_ODE_Solver_Peformance_for_Biological_Processes/sedml_import"



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
            #os.rename()
    except Exception as e:
        logger.warn(f"Failed to download sbml model from {sbml_url}, {e}.")
