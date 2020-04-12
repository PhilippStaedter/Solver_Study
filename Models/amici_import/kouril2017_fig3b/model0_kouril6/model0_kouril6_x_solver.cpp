#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kouril6(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ADP;
    x_solver[1] = ATP;
    x_solver[2] = BPG;
    x_solver[3] = P3G;
    x_solver[4] = pep;
    x_solver[5] = phos;
    x_solver[6] = pyr;
}