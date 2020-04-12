#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model1_saeidi1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = s1;
    x_solver[1] = s16;
    x_solver[2] = s17;
    x_solver[3] = s19;
    x_solver[4] = s2;
    x_solver[5] = s3;
    x_solver[6] = s4;
    x_solver[7] = s42;
    x_solver[8] = s45;
    x_solver[9] = s5;
}