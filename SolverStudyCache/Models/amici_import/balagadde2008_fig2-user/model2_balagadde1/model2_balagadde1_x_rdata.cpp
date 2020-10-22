#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model2_balagadde1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = A1;
    x_rdata[1] = A2;
    x_rdata[2] = C1;
    x_rdata[3] = C2;
    x_rdata[4] = IPTG;
    x_rdata[5] = sink;
    x_rdata[6] = source;
}