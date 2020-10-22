#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_montanez1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ARGex;
    x_rdata[1] = ARGin;
    x_rdata[2] = ORN;
}