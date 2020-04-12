#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_essunger3(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Ta;
    x_solver[1] = Tastarstar;
    x_solver[2] = Tm;
    x_solver[3] = Tmstarstar;
    x_solver[4] = Ttot;
    x_solver[5] = Tv;
    x_solver[6] = V;
}