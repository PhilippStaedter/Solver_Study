#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_piedrafita1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = S;
    x_solver[1] = ST;
    x_solver[2] = STU;
    x_solver[3] = STUS;
    x_solver[4] = STUST;
    x_solver[5] = STUSU;
    x_solver[6] = SU;
    x_solver[7] = SUST;
    x_solver[8] = SUSTU;
    x_solver[9] = T;
    x_solver[10] = U;
}