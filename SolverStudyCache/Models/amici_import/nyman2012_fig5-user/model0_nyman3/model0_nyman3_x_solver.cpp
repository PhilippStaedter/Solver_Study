#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_nyman3(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IR;
    x_solver[1] = IRS;
    x_solver[2] = IRSiP;
    x_solver[3] = IRi;
    x_solver[4] = IRiP;
    x_solver[5] = IRins;
    x_solver[6] = IRp;
    x_solver[7] = X;
    x_solver[8] = Xp;
}