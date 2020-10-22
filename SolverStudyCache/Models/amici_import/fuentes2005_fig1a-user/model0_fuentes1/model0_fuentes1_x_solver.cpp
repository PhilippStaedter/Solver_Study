#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_fuentes1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = e;
    x_solver[1] = ez;
    x_solver[2] = amici_w;
    x_solver[3] = z;
}