#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_esteban2013(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = M;
    x_solver[1] = MStar;
    x_solver[2] = T;
    x_solver[3] = TStar;
    x_solver[4] = V;
}