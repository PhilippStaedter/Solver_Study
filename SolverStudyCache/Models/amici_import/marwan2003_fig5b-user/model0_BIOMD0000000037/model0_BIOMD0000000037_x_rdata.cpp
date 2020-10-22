#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_BIOMD0000000037(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Gluc;
    x_rdata[1] = Pfr;
    x_rdata[2] = Pi;
    x_rdata[3] = Pr;
    x_rdata[4] = S;
    x_rdata[5] = V;
    x_rdata[6] = Xa;
    x_rdata[7] = Xi;
    x_rdata[8] = Ya;
    x_rdata[9] = Yi;
    x_rdata[10] = preS;
    x_rdata[11] = prepreS;
}