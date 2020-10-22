#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model2_balagadde1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = A1;
    x_solver[1] = A2;
    x_solver[2] = C1;
    x_solver[3] = C2;
    x_solver[4] = IPTG;
    x_solver[5] = sink;
    x_solver[6] = source;
}