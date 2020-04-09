# execute simulation script with settings of choice

######################################################################################
# very important: linSol, nonLinSolIter, solAlg must be their corresponding integers #
#                                                                                    #
# linSol:           Dense == 1  GMRES == 6  BiCGStab == 7   SPTFQMR == 8    KLU 00 9 #
# nonLinSolIter:    Functional == 1     Newton-type == 2                             #
# solAlg:           Adams == 1      BDF == 2                                         #
#                                                                                    #
######################################################################################

from SimAllSettings import *
from SimPriorSettings import *

# settings
atol = 1e-6
rtol = 1e-6
linSol = 9
nonLinSol = 2
solAlg = 2
maxStep = 10000


#simulate(atol, rtol, linSol, nonLinSol, solAlg, maxStep)
prior(atol, rtol, linSol, nonLinSol, solAlg, maxStep)
# useless = os.system('clear')