#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_fraser1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = a1;
    x_solver[1] = a2;
    x_solver[2] = a3;
    x_solver[3] = x1;
    x_solver[4] = x2;
    x_solver[5] = x3;
}