#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model3_chance3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = E;
    x_rdata[1] = P;
    x_rdata[2] = Q;
    x_rdata[3] = X;
}