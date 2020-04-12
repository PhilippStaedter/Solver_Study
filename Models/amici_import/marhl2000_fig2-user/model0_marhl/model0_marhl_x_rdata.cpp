#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_marhl(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CaER;
    x_rdata[1] = CaM;
    x_rdata[2] = CaPr;
    x_rdata[3] = Ca_cyt;
    x_rdata[4] = Pr;
}