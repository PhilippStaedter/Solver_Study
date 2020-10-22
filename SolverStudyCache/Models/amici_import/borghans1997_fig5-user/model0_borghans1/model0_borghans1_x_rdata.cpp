#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_borghans1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EC;
    x_rdata[1] = Fraction_Inactive_Channels;
    x_rdata[2] = Rho;
    x_rdata[3] = Y;
    x_rdata[4] = Z;
}