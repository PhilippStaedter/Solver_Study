# Solver_Study

This repository contains all python files, (supplementary) figures and references for the paper 'Assessment and Prediction of ODE Solver Performance for BiologicalProcesses'.

## Problems that can occur by downloading this repository on Windows:
	names for specific models can be too long to save (e.g. for the model 'Morris2002')
        Work-around ???

The study can be reproduced by following the order in which all scripts must be executed.
All python scripts can be found in 'Solver_Study/Python_Scripts'.
The study was done by using the software versions listed in 'Solver_Study/Software/requirements.txt'.

## 1 Create model collection 
### 1.1 Download all sedml and sbml models from the JWS Online Database
	script_download_all_sedml_models.py

### 1.2 Download chosen sbml models from the BioModels Database
	downloaded "by hand" --- no script available
	all those models are included in 'Solver_Study/Models/all_models' and their amici-imported-versions in 'Solver_Study/all_models/amici_import'

### 1.3 Import all sbml models to AMICI
	sbml2amici.py
	sbml2amici_BioModelsDatabase.py

### 1.4 Compare the state trajectories of the local simulation to the in-built simulation routine of JWS
	compareStateTrajectories_1.py
	compareStateTrajectories_2.py

### 1.5 Derive the whole model collection
	correctStateTrajectories.py

To skip step 1, the whole benchmark collection is available in 'Solver_Study/Models'

## 2 Solver settings study
### 2.1 Get all Data
	execute_WholeStudy.py
	execute_Tolerances.py

To skip step 2.1, all data files can be found in 'Solver_Study/Data'

### 2.2 Visualize all results according to the ordner seen in the paper
##### 2.2.1 Basic Properties
	    paper_plotCompareStateTrajectories.py
 	    paper_plotCompareStateTrajectories_2.py
	    paper_plotFirstResults.py
	    paper_plotFirstStudy.py

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


