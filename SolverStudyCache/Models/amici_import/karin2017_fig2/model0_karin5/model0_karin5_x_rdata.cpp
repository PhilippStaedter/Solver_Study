#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_karin5(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = beta;
    x_rdata[1] = g;
    x_rdata[2] = ins;
    x_rdata[3] = mbeta;
    x_rdata[4] = tamox;
}