#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_orfao1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = II;
    x_rdata[1] = IIa;
    x_rdata[2] = IIa_ATIII;
    x_rdata[3] = IIa_alpha2M;
    x_rdata[4] = PL;
    x_rdata[5] = PT;
    x_rdata[6] = RVV;
    x_rdata[7] = V;
    x_rdata[8] = Va;
    x_rdata[9] = X;
    x_rdata[10] = Xa;
    x_rdata[11] = Xa_ATIII;
}