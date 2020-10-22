#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_delcontezerail1_Fig4D(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = R5;
    x_solver[1] = R7;
    x_solver[2] = r5;
    x_solver[3] = r7;
}