#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_alexander1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = A;
    x_solver[1] = A_im;
    x_solver[2] = E;
    x_solver[3] = G;
    x_solver[4] = R;
}