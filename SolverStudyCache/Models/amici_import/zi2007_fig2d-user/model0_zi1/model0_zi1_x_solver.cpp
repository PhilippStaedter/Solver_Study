#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_zi1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = LRC_Cave;
    x_solver[1] = LRC_EE;
    x_solver[2] = LRC_Surf;
    x_solver[3] = Smad2c;
    x_solver[4] = Smad2n;
    x_solver[5] = Smad4c;
    x_solver[6] = Smad4n;
    x_solver[7] = Smads_Complex_c;
    x_solver[8] = Smads_Complex_n;
    x_solver[9] = T1R_Cave;
    x_solver[10] = T1R_EE;
    x_solver[11] = T1R_Surf;
    x_solver[12] = T2R_Cave;
    x_solver[13] = T2R_EE;
    x_solver[14] = T2R_Surf;
    x_solver[15] = TGF_beta;
}