#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model2_perelson3(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = T;
    x_solver[1] = Tstar;
    x_solver[2] = Tstarstar;
    x_solver[3] = Ttot;
    x_solver[4] = V;
}