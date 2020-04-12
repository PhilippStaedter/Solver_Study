#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_hockin2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = APC;
    x_rdata[1] = HC;
    x_rdata[2] = HC3;
    x_rdata[3] = HC36;
    x_rdata[4] = HC5;
    x_rdata[5] = HC53;
    x_rdata[6] = HC536;
    x_rdata[7] = HC56;
    x_rdata[8] = LC;
    x_rdata[9] = LC_APC;
    x_rdata[10] = Va;
    x_rdata[11] = Va3;
    x_rdata[12] = Va36;
    x_rdata[13] = Va36_APC;
    x_rdata[14] = Va3_APC;
    x_rdata[15] = Va5;
    x_rdata[16] = Va53;
    x_rdata[17] = Va536;
    x_rdata[18] = Va536_APC;
    x_rdata[19] = Va53_APC;
    x_rdata[20] = Va56;
    x_rdata[21] = Va56_APC;
    x_rdata[22] = Va5_APC;
    x_rdata[23] = VaA3;
    x_rdata[24] = VaA36;
    x_rdata[25] = VaA53;
    x_rdata[26] = VaA536;
    x_rdata[27] = VaLCA1;
    x_rdata[28] = VaLCA1_APC;
    x_rdata[29] = Va_APC;
}