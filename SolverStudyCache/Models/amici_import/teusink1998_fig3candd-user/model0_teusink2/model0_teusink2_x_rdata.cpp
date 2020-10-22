#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_teusink2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ATP;
    x_rdata[1] = Fru16P2;
    x_rdata[2] = Glc;
    x_rdata[3] = HMP;
}