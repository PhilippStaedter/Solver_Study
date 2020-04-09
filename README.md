# Solver_Study

This repository contains all python files, (supplementary) figures and references for the paper 'Assessment and Prediction of ODE Solver Performance for BiologicalProcesses'.

The study can be reproduced by following the order in which all scripts must be executed.

## 1 Create model collection 
### 1.1 Download all sedml and sbml models from the JWS Online Database
	script_download_all_sedml_models.py

### 1.2 Download chosen sbml models from the BioModels Database
	downloaded "by hand" --- no script available

### 1.3 Import all sbml models to AMICI
	sbml2amici.py
	sbml2amici_BioModelsDatabase.py

### 1.4 Check for correctly generated .so files and delete all other sbml models 
	no explicit script yet 

### 1.5 Compare the state trajectories of the local simulation to the in-built simulation routine of JWS
	compareStateTrajectories_1.py
	compareStateTrajectories_2.py

### 1.6 Derive the whole model collection
	correctStateTrajectories.py


## 2 Solver settings study
### 2.1 Get all Data
	execute_WholeStudy.py
	execute_Tolerances.py

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

##### 2.2.5 Multistep method study
	    paper_plotAdamsBDF.py
	    paper_plotAdamsBDF_2.py	



