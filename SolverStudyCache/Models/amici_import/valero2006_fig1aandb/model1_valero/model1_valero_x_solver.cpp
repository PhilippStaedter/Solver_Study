#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_valero(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ADP;
    x_solver[1] = AMP;
    x_solver[2] = ATP;
    x_solver[3] = Lac;
    x_solver[4] = NADH;
    x_solver[5] = Pyr;
}