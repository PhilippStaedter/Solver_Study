#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_bindschadler1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = c1;
    x_rdata[1] = c2;
    x_rdata[2] = h1;
    x_rdata[3] = h2;
}