#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_tripathi1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = A;
    x_rdata[1] = EXT;
    x_rdata[2] = I1;
    x_rdata[3] = I2;
    x_rdata[4] = S;
}