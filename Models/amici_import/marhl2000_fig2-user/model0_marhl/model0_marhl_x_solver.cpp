#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_marhl(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = CaER;
    x_solver[1] = CaM;
    x_solver[2] = CaPr;
    x_solver[3] = Ca_cyt;
    x_solver[4] = Pr;
}