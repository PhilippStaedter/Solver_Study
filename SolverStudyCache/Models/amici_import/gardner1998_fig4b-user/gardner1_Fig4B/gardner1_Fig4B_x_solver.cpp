#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_gardner1_Fig4B(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = C;
    x_solver[1] = M;
    x_solver[2] = X;
    x_solver[3] = Y;
    x_solver[4] = Z;
}