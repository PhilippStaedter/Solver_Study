#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_tiago1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = s1;
    x_solver[1] = s10;
    x_solver[2] = s11;
    x_solver[3] = s12;
    x_solver[4] = s13;
    x_solver[5] = s14;
    x_solver[6] = s15;
    x_solver[7] = s16;
    x_solver[8] = s17;
    x_solver[9] = s2;
    x_solver[10] = s3;
    x_solver[11] = s4;
    x_solver[12] = s5;
    x_solver[13] = s6;
    x_solver[14] = s7;
    x_solver[15] = s8;
    x_solver[16] = s9;
}