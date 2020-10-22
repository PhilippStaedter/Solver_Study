#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_bier2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = G1;
    x_rdata[1] = G2;
    x_rdata[2] = T1;
    x_rdata[3] = T2;
}