#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_BIOMD0000000037(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Gluc;
    x_solver[1] = Pfr;
    x_solver[2] = Pi;
    x_solver[3] = Pr;
    x_solver[4] = S;
    x_solver[5] = V;
    x_solver[6] = Xa;
    x_solver[7] = Xi;
    x_solver[8] = Ya;
    x_solver[9] = Yi;
    x_solver[10] = preS;
    x_solver[11] = prepreS;
}