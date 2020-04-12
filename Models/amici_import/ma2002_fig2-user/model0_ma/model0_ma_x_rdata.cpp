#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_ma(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ACA;
    x_rdata[1] = CAR1;
    x_rdata[2] = ERK2;
    x_rdata[3] = PKA;
    x_rdata[4] = REGA;
    x_rdata[5] = excAMP;
    x_rdata[6] = incAMP;
}