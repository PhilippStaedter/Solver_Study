#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_miao2008(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = T;
    x_rdata[1] = Tm;
    x_rdata[2] = Tmw;
    x_rdata[3] = Tw;
}