#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_lee2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = E;
    x_rdata[1] = E_M;
    x_rdata[2] = E_P2;
    x_rdata[3] = E_P_1;
    x_rdata[4] = E_P_2;
    x_rdata[5] = M;
    x_rdata[6] = P;
    x_rdata[7] = P2;
    x_rdata[8] = T;
}