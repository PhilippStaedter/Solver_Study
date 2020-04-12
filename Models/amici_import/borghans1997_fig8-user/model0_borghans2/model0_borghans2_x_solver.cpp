#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_borghans2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = A;
    x_solver[1] = EC;
    x_solver[2] = Y;
    x_solver[3] = Z;
}