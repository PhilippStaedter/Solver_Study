#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model22_karin1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = amici_y;
    x_rdata[1] = z;
}