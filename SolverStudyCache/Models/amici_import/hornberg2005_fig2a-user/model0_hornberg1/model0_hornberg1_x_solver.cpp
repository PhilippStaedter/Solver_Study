#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_hornberg1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Inh;
    x_solver[1] = R;
    x_solver[2] = Rin;
    x_solver[3] = x1;
    x_solver[4] = x1p;
    x_solver[5] = x2;
    x_solver[6] = x2p;
    x_solver[7] = x3;
    x_solver[8] = x3p;
}