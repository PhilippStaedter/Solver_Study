# download all models from the JWS database

from sedml_import import *
from logging_util import *
import logging


log_to_console(level=logging.INFO)
log_to_file(level=logging.WARN, filename="../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/import_warnings.log")
logger = logging.getLogger(LOGGER_NAME)
logger.setLevel(logging.INFO)

# download all models
download_all_sedml_models_from_jws()
