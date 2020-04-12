#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_xie1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CC;
    x_rdata[1] = CCPT;
    x_rdata[2] = CLK;
    x_rdata[3] = CYC;
    x_rdata[4] = PDP;
    x_rdata[5] = PER;
    x_rdata[6] = PT;
    x_rdata[7] = TIM;
    x_rdata[8] = VRI;
    x_rdata[9] = clkm;
    x_rdata[10] = clkp;
    x_rdata[11] = pdpm;
    x_rdata[12] = pdpp;
    x_rdata[13] = perm;
    x_rdata[14] = perp;
    x_rdata[15] = prcpdp;
    x_rdata[16] = prcper;
    x_rdata[17] = prct;
    x_rdata[18] = prcv;
    x_rdata[19] = prpc;
    x_rdata[20] = prvc;
    x_rdata[21] = timm;
    x_rdata[22] = timp;
    x_rdata[23] = vrim;
    x_rdata[24] = vrip;
}