#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_dupont2_Fig4B(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = W_star;
    x_rdata[1] = Wt;
    x_rdata[2] = Y;
    x_rdata[3] = Z;
}