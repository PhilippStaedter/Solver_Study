#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_bier2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = G1;
    x_solver[1] = G2;
    x_solver[2] = T1;
    x_solver[3] = T2;
}