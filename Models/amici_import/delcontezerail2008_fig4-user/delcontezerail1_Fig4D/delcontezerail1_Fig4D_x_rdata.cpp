#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_delcontezerail1_Fig4D(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = R5;
    x_rdata[1] = R7;
    x_rdata[2] = r5;
    x_rdata[3] = r7;
}