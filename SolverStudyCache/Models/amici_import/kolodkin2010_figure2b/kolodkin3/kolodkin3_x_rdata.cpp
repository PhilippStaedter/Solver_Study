#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_kolodkin3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Lc;
    x_rdata[1] = Ln_;
    x_rdata[2] = NRLc;
    x_rdata[3] = NRLm;
    x_rdata[4] = NRLn;
    x_rdata[5] = NRc;
    x_rdata[6] = NRm;
    x_rdata[7] = NRn;
    x_rdata[8] = RE;
    x_rdata[9] = REL;
}