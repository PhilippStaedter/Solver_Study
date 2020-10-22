#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_lee4(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = E;
    x_rdata[1] = E_M;
    x_rdata[2] = E_M1;
    x_rdata[3] = E_P1;
    x_rdata[4] = E_P2;
    x_rdata[5] = E_P21;
    x_rdata[6] = E_P_1;
    x_rdata[7] = E_P_2;
    x_rdata[8] = M;
    x_rdata[9] = M1;
    x_rdata[10] = P;
    x_rdata[11] = P2;
    x_rdata[12] = P21;
    x_rdata[13] = T;
}