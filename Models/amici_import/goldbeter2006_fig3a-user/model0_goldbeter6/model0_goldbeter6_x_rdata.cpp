#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_goldbeter6(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = P;
    x_rdata[1] = Q;
    x_rdata[2] = R;
}