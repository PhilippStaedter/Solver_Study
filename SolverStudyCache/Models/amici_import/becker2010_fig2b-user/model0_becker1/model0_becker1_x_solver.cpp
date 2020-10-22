#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_becker1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Epo;
    x_solver[1] = EpoR;
    x_solver[2] = Epo_EpoR;
    x_solver[3] = Epo_EpoRi;
    x_solver[4] = dEpoe;
    x_solver[5] = dEpoi;
}