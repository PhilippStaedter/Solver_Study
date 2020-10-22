#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_aguda1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CYCDcdk4;
    x_rdata[1] = CYCDcdk4p27;
    x_rdata[2] = CYCEcdk2p27;
    x_rdata[3] = E2F;
    x_rdata[4] = aCYCEcdk2;
    x_rdata[5] = aCYCEcdk20;
    x_rdata[6] = aCYCEcdk21;
    x_rdata[7] = iCYCEcdk2;
    x_rdata[8] = p27;
    x_rdata[9] = pRB;
    x_rdata[10] = pRBE2F;
    x_rdata[11] = xvar;
}