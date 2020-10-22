#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kouril8(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ADP;
    x_rdata[1] = ATP;
    x_rdata[2] = BPG;
    x_rdata[3] = GA;
    x_rdata[4] = P3G;
    x_rdata[5] = gap;
    x_rdata[6] = glc;
    x_rdata[7] = nadp;
    x_rdata[8] = nadph;
    x_rdata[9] = pep;
    x_rdata[10] = pi;
    x_rdata[11] = pyr;
    x_rdata[12] = sink;
}