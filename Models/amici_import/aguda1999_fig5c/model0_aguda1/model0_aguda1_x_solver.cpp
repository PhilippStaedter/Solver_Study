#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_aguda1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CYCDcdk4;
    x_solver[1] = CYCDcdk4p27;
    x_solver[2] = CYCEcdk2p27;
    x_solver[3] = E2F;
    x_solver[4] = aCYCEcdk2;
    x_solver[5] = aCYCEcdk20;
    x_solver[6] = aCYCEcdk21;
    x_solver[7] = iCYCEcdk2;
    x_solver[8] = p27;
    x_solver[9] = pRB;
    x_solver[10] = pRBE2F;
    x_solver[11] = xvar;
}