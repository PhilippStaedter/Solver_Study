#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_teusink2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ATP;
    x_solver[1] = Fru16P2;
    x_solver[2] = Glc;
    x_solver[3] = HMP;
}