#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_fuentes1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = e;
    x_rdata[1] = ez;
    x_rdata[2] = amici_w;
    x_rdata[3] = z;
}