#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model3_chance3(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = E;
    x_solver[1] = P;
    x_solver[2] = Q;
    x_solver[3] = X;
}