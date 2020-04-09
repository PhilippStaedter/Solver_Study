# SBML2AMICI
import libsbml
import importlib
import amici
import os
import sys
import numpy as np
import logging
import shutil
import pandas as pd


# create .tsv file
tsv_table = pd.DataFrame(columns=['id', 'states', 'reactions', 'error_message'])      # index=range() can be left out

# important paths
models_path = "../sbml2amici/amici_models"
models_base_path = "../sbml2amici"
base_path = "./BioModelsDatabase_models"


# create directory for all future amici models
if not os.path.exists(models_path):
    os.makedirs(models_path)

# create logger object
logger = logging.getLogger()

# initialize the log settings
logging.basicConfig(filename='all_logs',level=logging.DEBUG)

# list of all directories + SBML files
list_directory = os.listdir(base_path)
list_directory = sorted(list_directory)

# set row-counter for .tsv file
counter = 0

for models in list_directory:

    models = 'Bungay2003'

    list_files = os.listdir(base_path + '/' + models + '/sbml_models')
    list_files = sorted(list_files)

    for files in list_files:
        sbml_file = base_path + '/' + models + '/sbml_models/' + files
        model_name, other_stuff = files.split(".",1)
        model_output_dir = models_path + '/' + models + '/' + model_name

        # get new_observables()

        try:
            ## entry = {'id': None, 'states': None, 'reactions': None, 'parameters': None, 'error_message': None}
            # Append additional row in .tsv file
            tsv_table = tsv_table.append({}, ignore_index=True)

            # name
            tsv_table.loc[counter].id = '{' + models + '}' + '_' + '{' + files + '}'

            # read accompanying sbml file
            file = libsbml.readSBML(sbml_file)
            all_properties = file.getModel()
            num_states = len(all_properties.getListOfSpecies())
            num_reactions = len(all_properties.getListOfReactions())

            # Fill in .tsv file
            tsv_table.loc[counter].states = num_states
            tsv_table.loc[counter].reactions = num_reactions

            # Create SBML importer
            sbml_importer = amici.SbmlImporter(sbml_file)

            # SBML2AMICI
            sbml_importer.sbml2amici(model_name,
                        model_output_dir,
                        verbose=False)                                                                  # TRUE instead of FALSE?

            # read data out
            # sbml_importer.

            # Write 'OK' in 'error_message' coloumn
            tsv_table.loc[counter].error_message = 'OK'

            # increase counter by 1
            counter = counter + 1

        except Exception as e:
            error_info = str(e)
            # error_info = sys.exc_info()[0]
            print(error_info)
            logging.exception('Model failed: %s, %s', models, files)
            logging.info('\n')

            # Write the error message in 'error_message' coloumn
            tsv_table.loc[counter].error_message = error_info
            ## tsv_table.append(entry)
            # increase counter by 1
            counter = counter + 1

            # logging.exception(str(Exception))
            # except Exception as e:
            # logging.exeption(str(e), exc_info=True)
            # continue

# print tsv_table + save it
# print(tsv_table)
tsv_table.to_csv(path_or_buf=models_base_path + '/table_BioModelsDatabase.tsv', sep='\t', index=True)

# copy file 'all_logs' in new directory 'sbml2amici'
old_path = './all_logs'
new_path = models_base_path + '/all_logs'
shutil.move(old_path, new_path)