#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kouril6(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ADP;
    x_rdata[1] = ATP;
    x_rdata[2] = BPG;
    x_rdata[3] = P3G;
    x_rdata[4] = pep;
    x_rdata[5] = phos;
    x_rdata[6] = pyr;
}