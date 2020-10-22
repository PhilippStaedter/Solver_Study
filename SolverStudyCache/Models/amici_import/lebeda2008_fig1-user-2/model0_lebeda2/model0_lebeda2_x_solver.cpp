#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_lebeda2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = bound;
    x_solver[1] = free;
    x_solver[2] = lytic;
    x_solver[3] = translocate;
}