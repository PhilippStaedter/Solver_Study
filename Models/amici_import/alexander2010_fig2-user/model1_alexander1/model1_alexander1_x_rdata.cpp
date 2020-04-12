#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_alexander1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = A;
    x_rdata[1] = A_im;
    x_rdata[2] = E;
    x_rdata[3] = G;
    x_rdata[4] = R;
}