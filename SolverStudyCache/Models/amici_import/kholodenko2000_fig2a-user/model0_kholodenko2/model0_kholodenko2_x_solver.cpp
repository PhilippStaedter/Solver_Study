#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_kholodenko2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = MAPK;
    x_solver[1] = MAPK_P;
    x_solver[2] = MAPK_PP;
    x_solver[3] = MKK;
    x_solver[4] = MKKK;
    x_solver[5] = MKKK_P;
    x_solver[6] = MKK_P;
    x_solver[7] = MKK_PP;
}