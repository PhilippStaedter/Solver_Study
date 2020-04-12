#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_lee3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = II;
    x_rdata[1] = IIa;
    x_rdata[2] = M;
    x_rdata[3] = P2;
}