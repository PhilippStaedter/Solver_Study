#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model9_panteleev3(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = TFPI;
    x_solver[1] = VIIa_TF;
    x_solver[2] = VIIa_TF_X;
    x_solver[3] = VIIa_TF_Xa;
    x_solver[4] = X;
    x_solver[5] = Xa;
    x_solver[6] = Xa_TFPI;
    x_solver[7] = Xa_TFPI_VIIa_TF;
}