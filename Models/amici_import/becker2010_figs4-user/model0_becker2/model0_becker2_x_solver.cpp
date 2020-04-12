#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_becker2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EpoR;
    x_solver[1] = SAv;
    x_solver[2] = SAv_EpoR;
    x_solver[3] = SAv_EpoRi;
    x_solver[4] = dSAve;
    x_solver[5] = dSAvi;
}