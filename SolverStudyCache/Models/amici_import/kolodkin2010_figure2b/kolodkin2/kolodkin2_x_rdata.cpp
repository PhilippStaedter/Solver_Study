#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_kolodkin2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Lc;
    x_rdata[1] = Ln_;
    x_rdata[2] = NRLn;
    x_rdata[3] = NRn;
    x_rdata[4] = RE;
    x_rdata[5] = REL;
}