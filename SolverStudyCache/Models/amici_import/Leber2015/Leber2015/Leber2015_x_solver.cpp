#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Leber2015(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Cdiff;
    x_solver[1] = Commensal_Beneficial;
    x_solver[2] = Commensal_Dead;
    x_solver[3] = tDC_LP;
    x_solver[4] = tDC_MLN;
    x_solver[5] = Commensal_Harmful;
    x_solver[6] = N_Lum;
    x_solver[7] = E;
    x_solver[8] = E_d;
    x_solver[9] = iDC_E;
    x_solver[10] = E_i;
    x_solver[11] = M_LP;
    x_solver[12] = eDC_LP;
    x_solver[13] = M0;
    x_solver[14] = N_LP;
    x_solver[15] = Th17_LP;
    x_solver[16] = Th1_LP;
    x_solver[17] = iTreg_LP;
    x_solver[18] = eDC_MLN;
    x_solver[19] = iTreg_MLN;
    x_solver[20] = nT;
    x_solver[21] = Th17_MLN;
    x_solver[22] = Th1_MLN;
}