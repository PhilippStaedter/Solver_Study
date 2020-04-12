#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_zhao1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EXT;
    x_rdata[1] = I;
    x_rdata[2] = Ip;
    x_rdata[3] = S;
    x_rdata[4] = Sp;
}