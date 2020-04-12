#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_gardner1_Fig4B(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = C;
    x_rdata[1] = M;
    x_rdata[2] = X;
    x_rdata[3] = Y;
    x_rdata[4] = Z;
}