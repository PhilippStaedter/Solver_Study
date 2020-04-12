#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_hornberg1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Inh;
    x_rdata[1] = R;
    x_rdata[2] = Rin;
    x_rdata[3] = x1;
    x_rdata[4] = x1p;
    x_rdata[5] = x2;
    x_rdata[6] = x2p;
    x_rdata[7] = x3;
    x_rdata[8] = x3p;
}