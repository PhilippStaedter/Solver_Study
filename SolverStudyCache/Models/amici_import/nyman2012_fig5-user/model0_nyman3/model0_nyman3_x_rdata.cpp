#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_nyman3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = IR;
    x_rdata[1] = IRS;
    x_rdata[2] = IRSiP;
    x_rdata[3] = IRi;
    x_rdata[4] = IRiP;
    x_rdata[5] = IRins;
    x_rdata[6] = IRp;
    x_rdata[7] = X;
    x_rdata[8] = Xp;
}