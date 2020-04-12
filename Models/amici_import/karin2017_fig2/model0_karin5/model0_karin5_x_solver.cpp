#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_karin5(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = beta;
    x_solver[1] = g;
    x_solver[2] = ins;
    x_solver[3] = mbeta;
    x_solver[4] = tamox;
}