#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_kolodkin2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Lc;
    x_solver[1] = Ln_;
    x_solver[2] = NRLn;
    x_solver[3] = NRn;
    x_solver[4] = RE;
    x_solver[5] = REL;
}