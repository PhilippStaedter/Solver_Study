#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_orfao1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = II;
    x_solver[1] = IIa;
    x_solver[2] = IIa_ATIII;
    x_solver[3] = IIa_alpha2M;
    x_solver[4] = PL;
    x_solver[5] = PT;
    x_solver[6] = RVV;
    x_solver[7] = V;
    x_solver[8] = Va;
    x_solver[9] = X;
    x_solver[10] = Xa;
    x_solver[11] = Xa_ATIII;
}