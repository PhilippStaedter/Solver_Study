#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_essunger1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Ta;
    x_rdata[1] = Tastarstar;
    x_rdata[2] = Tm;
    x_rdata[3] = Tmstarstar;
    x_rdata[4] = Ttot;
    x_rdata[5] = Tv;
    x_rdata[6] = V;
}