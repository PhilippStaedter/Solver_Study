#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_lou1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Ib;
    x_rdata[1] = It;
    x_rdata[2] = Iv;
    x_rdata[3] = Sb;
    x_rdata[4] = St;
    x_rdata[5] = Sv;
}