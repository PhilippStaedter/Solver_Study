#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_lee4(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = E;
    x_solver[1] = E_M;
    x_solver[2] = E_M1;
    x_solver[3] = E_P1;
    x_solver[4] = E_P2;
    x_solver[5] = E_P21;
    x_solver[6] = E_P_1;
    x_solver[7] = E_P_2;
    x_solver[8] = M;
    x_solver[9] = M1;
    x_solver[10] = P;
    x_solver[11] = P2;
    x_solver[12] = P21;
    x_solver[13] = T;
}