#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_fisher1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Act_C_Cyt;
    x_solver[1] = Act_C_Nuc;
    x_solver[2] = Ca_Cyt;
    x_solver[3] = Ca_Nuc;
    x_solver[4] = Inact_C_Cyt;
    x_solver[5] = Inact_C_Nuc;
    x_solver[6] = NFAT_Act_C_Cyt;
    x_solver[7] = NFAT_Act_C_Nuc;
    x_solver[8] = NFAT_Cyt;
    x_solver[9] = NFAT_Nuc;
    x_solver[10] = NFAT_Pi_Act_C_Cyt;
    x_solver[11] = NFAT_Pi_Act_C_Nuc;
    x_solver[12] = NFAT_Pi_Cyt;
    x_solver[13] = NFAT_Pi_Nuc;
}