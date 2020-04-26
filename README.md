# Solver_Study

This repository contains all python files, (supplementary) figures and references for the paper 'Assessment and Prediction of ODE Solver Performance for BiologicalProcesses'.

## 0 Problems that can oocur
### 0.1 by downloading this repository on Windows:
	names for specific models can be too long to save (error could occur e.g. for the model 'Morris2002')
        Work-around: Omit this/these model(s)

### 0.2 by working with python packages
	if the package tqdm does not show a progress bar, the button "Emulate terminal in output console" in pythons' configurations must be siwtched on.

### 0.3 by doing step 1 only partially
	in step 1, the main folder 'Assessment_of_ODE_Solver_Performance_for_Biological_Processes' with its subfolders 'BioModelsDatabase_models', 'json_files', 'sbml2amici' and 'sedml_models' will be 	 created. If some but not all steps in 1 are executed, the following scripts in step 2 might not work as intended. To some extend, all scripts cover many possibilities but maybe not all.
	Hence, to guarantee a successful reproduction of the study, please execute step 1 with all substeps or not at all.

The study can be reproduced by following the order in which all scripts have to be executed.
All python scripts can be found in 'Solver_Study/Python_Scripts'.
The study was done by using the software versions listed in 'Solver_Study/Software/requirements.txt'.
By reproducing this study, all results are stored in a therefore created folder 'Assessment_of_ODE_Solver_Performance_for_Biological_Processes' on the same level as the downloaded repository.

## 1 Create model collection 
### 1.1 Download all sedml and sbml models from the JWS Online Database
	script_download_all_sedml_models.py

### 1.2 Download chosen sbml models from the BioModels Database
	script_download_all_sbml_biomodels.py

### 1.3 Import all sbml models to AMICI
	sbml2amici.py
	sbml2amici_BioModelsDatabase.py

### 1.4 Compare the state trajectories of the local simulation to the in-built simulation routine of JWS
	compareStateTrajectories_1.py
	compareStateTrajectories_2.py

### 1.5 Derive the whole model collection
	correctStateTrajectories.py

To skip step 1, the whole benchmark collection is available in 'Solver_Study/Models'.
If step 1 was skipped, the main folder 'Assessment_of_ODE_Solver_Performance_for_Biological_Processes' or the subfolders 'BioModelsDatabase_models', 'json_files', 'sbml2amici' and 'sedml_models' will not exist. 
In this case, the next functions will automatically take the results from 'Solver_Study/Models' of the repository.
Additionally, the main folder 'Assessment_of_ODE_Solver_Performance_for_Biological_Processes' will be created. 

## 2 Solver settings study
### 2.1 Get all Data
	execute_WholeStudy.py
	execute_Tolerances.py

To skip step 2.1, all data files can be found in 'Solver_Study/Data'
If step 2.1 was skipped, the subfolder 'Data' will not exist.
If step 1 was skipped, the subfolder 'json_files' will not exist.
In this case, the next functions will automatically take the results from 'Solver_Study/Data' of the repository.
Additionally, the main folder 'Assessment_of_ODE_Solver_Performance_for_Biological_Processes' will be created. 

### 2.2 Visualize all results according to the ordner seen in the paper
Important: All figures created by the following scripts are not stored automatically!
##### 2.2.1 Basic Properties
	    paper_plotCompareStateTrajectories.py (not needed)
 	    paper_plotCompareStateTrajectories_2.py (Supplementary)
	    paper_plotFirstResults.py (not needed)
	    paper_plotFirstStudy.py (Main Manuskript)

##### 2.2.2 Non-linear solver study
	    paper_plotNonLinSol.py
	    paper_plotNonLinSol_2.py
	    paper_plotNonLinSol_3.py
	    paper_plotNonLinSol_4.py

##### 2.2.3 Linear solver study
	    paper_plotLinearSolver.py
	    paper_plotScatter.py
	    paper_plotScatter_2.py
	    paper_plotScatter_3.py
	    paper_plotScatter_4.py	

##### 2.2.4 Tolerances study
	    paper_plotBoxPlot.py
	    paper_plotBoxPlot_2.py
	    paper_plotHistogram.py

##### 2.2.5 Solver algorithm study
	    paper_plotAdamsBDF.py
	    paper_plotAdamsBDF_2.py	

To skip step 2.2, all figures can be found in 'Solver_Study/Figures'


