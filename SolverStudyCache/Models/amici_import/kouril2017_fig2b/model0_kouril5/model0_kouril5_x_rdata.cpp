#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kouril5(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = BPG;
    x_rdata[1] = P3G;
}