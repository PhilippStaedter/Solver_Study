#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Yang2007(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = x1;
    x_solver[1] = x10;
    x_solver[2] = x11;
    x_solver[3] = x12;
    x_solver[4] = x13;
    x_solver[5] = x14;
    x_solver[6] = x15;
    x_solver[7] = x16;
    x_solver[8] = x17;
    x_solver[9] = x18;
    x_solver[10] = x19;
    x_solver[11] = x2;
    x_solver[12] = x20;
    x_solver[13] = x21;
    x_solver[14] = x22;
    x_solver[15] = x23;
    x_solver[16] = x24;
    x_solver[17] = x25;
    x_solver[18] = x3;
    x_solver[19] = x4;
    x_solver[20] = x5;
    x_solver[21] = x6;
    x_solver[22] = x7;
    x_solver[23] = x8;
    x_solver[24] = x9;
}