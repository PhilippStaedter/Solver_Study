#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_fraser1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = a1;
    x_rdata[1] = a2;
    x_rdata[2] = a3;
    x_rdata[3] = x1;
    x_rdata[4] = x2;
    x_rdata[5] = x3;
}