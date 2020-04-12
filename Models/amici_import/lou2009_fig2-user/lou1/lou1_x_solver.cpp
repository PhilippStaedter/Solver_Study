#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_lou1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Ib;
    x_solver[1] = It;
    x_solver[2] = Iv;
    x_solver[3] = Sb;
    x_solver[4] = St;
    x_solver[5] = Sv;
}