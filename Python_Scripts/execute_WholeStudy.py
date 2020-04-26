# execute simulation script with all solver setiing combinations

######################################################################################
# very important: linSol, nonLinSolIter, solAlg must be their corresponding integers #
#                                                                                    #
# linSol:           Dense == 1  GMRES == 6  BiCGStab == 7   SPTFQMR == 8    KLU == 9 #
# nonLinSolIter:    Functional == 1     Newton-type == 2                             #
# solAlg:           Adams == 1      BDF == 2                                         #
#                                                                                    #
######################################################################################

from SimAllSettings import *
from tqdm import trange
import os


# settings
atol = [1e-8, 1e-6, 1e-12, 1e-10, 1e-14, 1e-16, 1e-8]
rtol = [1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-8, 1e-16]
linSol = [1, 6, 7, 8, 9]
iter = [1, 2]
solAlg = [1, 2]
maxStep = 10000
study_indicator = 1                                         # for the main study, don't change!
skip_indicator = 0                                          # tracks whether the step 1 was skipped (0) or not (1)
if os.path.exists('../../Assessment_of_ODE_Solver_Performance_for_Biological_Processes/sedml_models'):
    skip_indicator = 1

# run simulation
for iSolALg in trange(0, len(solAlg)):
    for iNonLin in trange(0, len(iter)):
        for iLinSol in trange(0, len(linSol)):
            for iToleranceTupel in trange(0, len(atol)):
                simulate(atol[iToleranceTupel], rtol[iToleranceTupel], linSol[iLinSol], iter[iNonLin], solAlg[iSolALg], maxStep, study_indicator, skip_indicator)
                useless_variable = os.system('clear')
                pass
