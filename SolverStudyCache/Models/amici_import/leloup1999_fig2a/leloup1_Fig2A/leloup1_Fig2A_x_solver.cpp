#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_leloup1_Fig2A(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CC;
    x_solver[1] = Cn;
    x_solver[2] = Mp;
    x_solver[3] = Mt;
    x_solver[4] = P0;
    x_solver[5] = P1;
    x_solver[6] = P2;
    x_solver[7] = T0;
    x_solver[8] = T1;
    x_solver[9] = T2;
}