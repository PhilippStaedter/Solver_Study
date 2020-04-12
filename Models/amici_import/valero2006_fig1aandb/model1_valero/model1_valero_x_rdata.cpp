#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_valero(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ADP;
    x_rdata[1] = AMP;
    x_rdata[2] = ATP;
    x_rdata[3] = Lac;
    x_rdata[4] = NADH;
    x_rdata[5] = Pyr;
}