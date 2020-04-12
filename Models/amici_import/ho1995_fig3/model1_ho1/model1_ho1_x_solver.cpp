#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_ho1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = T;
}