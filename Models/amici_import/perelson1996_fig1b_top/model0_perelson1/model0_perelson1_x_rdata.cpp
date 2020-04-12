#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_perelson1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Tstar;
    x_rdata[1] = V;
    x_rdata[2] = Vin;
    x_rdata[3] = Vni;
}