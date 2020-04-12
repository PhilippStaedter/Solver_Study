#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kouril8(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ADP;
    x_solver[1] = ATP;
    x_solver[2] = BPG;
    x_solver[3] = GA;
    x_solver[4] = P3G;
    x_solver[5] = gap;
    x_solver[6] = glc;
    x_solver[7] = nadp;
    x_solver[8] = nadph;
    x_solver[9] = pep;
    x_solver[10] = pi;
    x_solver[11] = pyr;
    x_solver[12] = sink;
}