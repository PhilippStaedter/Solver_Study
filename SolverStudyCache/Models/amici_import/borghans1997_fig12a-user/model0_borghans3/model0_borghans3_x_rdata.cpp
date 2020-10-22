#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_borghans3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EC;
    x_rdata[1] = X;
    x_rdata[2] = Y;
    x_rdata[3] = Z;
}