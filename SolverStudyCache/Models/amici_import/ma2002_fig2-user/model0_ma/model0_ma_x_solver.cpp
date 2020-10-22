#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_ma(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ACA;
    x_solver[1] = CAR1;
    x_solver[2] = ERK2;
    x_solver[3] = PKA;
    x_solver[4] = REGA;
    x_solver[5] = excAMP;
    x_solver[6] = incAMP;
}