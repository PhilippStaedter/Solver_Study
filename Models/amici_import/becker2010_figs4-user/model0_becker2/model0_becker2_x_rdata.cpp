#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_becker2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EpoR;
    x_rdata[1] = SAv;
    x_rdata[2] = SAv_EpoR;
    x_rdata[3] = SAv_EpoRi;
    x_rdata[4] = dSAve;
    x_rdata[5] = dSAvi;
}