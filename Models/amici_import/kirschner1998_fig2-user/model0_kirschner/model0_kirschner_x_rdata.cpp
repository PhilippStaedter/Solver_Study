#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kirschner(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = T;
    x_rdata[1] = V;
}