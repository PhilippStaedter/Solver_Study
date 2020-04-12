#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_kolodkin5(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Lc;
    x_rdata[1] = Ln_;
    x_rdata[2] = NR;
    x_rdata[3] = NRLc;
    x_rdata[4] = NRLn;
    x_rdata[5] = NRc;
    x_rdata[6] = RE;
    x_rdata[7] = REL;
}