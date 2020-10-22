#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_hockin2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*Va*k1;
    dwdx[1] = 1.0*Va3*k1;
    dwdx[2] = 1.0*LC*k1;
    dwdx[3] = 1.0*Va5*k1;
    dwdx[4] = 1.0*VaLCA1*k1;
    dwdx[5] = 1.0*Va53*k1;
    dwdx[6] = 1.0*Va56*k1;
    dwdx[7] = 1.0*Va36*k1;
    dwdx[8] = 1.0*Va536*k1;
    dwdx[9] = -1.0*LC*k10;
    dwdx[10] = -1.0*LC*k10;
    dwdx[11] = -1.0*LC*k10;
    dwdx[12] = -1.0*LC*k10;
    dwdx[13] = -1.0*LC*k10;
    dwdx[14] = -1.0*LC*k10;
    dwdx[15] = -1.0*LC*k10;
    dwdx[16] = -1.0*HC*k10;
    dwdx[17] = -1.0*HC5*k10;
    dwdx[18] = -1.0*HC3*k10;
    dwdx[19] = -1.0*HC56*k10;
    dwdx[20] = -1.0*HC53*k10;
    dwdx[21] = -1.0*HC36*k10;
    dwdx[22] = -1.0*HC536*k10;
    dwdx[23] = 1.0*APC*k1;
    dwdx[24] = -1.0*k2;
    dwdx[25] = 1.0*APC*k1;
    dwdx[26] = 1.0*k9;
    dwdx[27] = 1.0*k9;
    dwdx[28] = 1.0*APC*k1;
    dwdx[29] = 1.0*k7;
    dwdx[30] = 1.0*k9;
    dwdx[31] = 1.0*k7;
    dwdx[32] = 1.0*APC*k1;
    dwdx[33] = 1.0*k5;
    dwdx[34] = 1.0*k7;
    dwdx[35] = -1.0*k2;
    dwdx[36] = 1.0*k6;
    dwdx[37] = 1.0*k5;
    dwdx[38] = -1.0*k2;
    dwdx[39] = 1.0*k7;
    dwdx[40] = 1.0*k9;
    dwdx[41] = 1.0*APC*k1;
    dwdx[42] = 1.0*k9;
    dwdx[43] = 1.0*k7;
    dwdx[44] = 1.0*APC*k1;
    dwdx[45] = 1.0*k9;
    dwdx[46] = 1.0*k7;
    dwdx[47] = 1.0*APC*k1;
    dwdx[48] = 1.0*k7;
    dwdx[49] = -1.0*k2;
    dwdx[50] = 1.0*k6;
    dwdx[51] = 1.0*k7;
    dwdx[52] = -1.0*k2;
    dwdx[53] = 1.0*k9;
    dwdx[54] = 1.0*APC*k1;
    dwdx[55] = 1.0*k3;
    dwdx[56] = -1.0*k2;
    dwdx[57] = 1.0*k3;
    dwdx[58] = -1.0*k2;
    dwdx[59] = 1.0*k6;
    dwdx[60] = -1.0*VaLCA1*k8;
    dwdx[61] = -1.0*VaLCA1_APC*k8;
    dwdx[62] = -1.0*VaLCA1*k8;
    dwdx[63] = -1.0*VaLCA1_APC*k8;
    dwdx[64] = -1.0*VaLCA1*k8;
    dwdx[65] = -1.0*VaLCA1_APC*k8;
    dwdx[66] = -1.0*VaLCA1*k8;
    dwdx[67] = -1.0*VaLCA1_APC*k7;
    dwdx[68] = -1.0*VaA3*k8;
    dwdx[69] = -1.0*VaA53*k8;
    dwdx[70] = -1.0*VaA36*k8;
    dwdx[71] = -1.0*VaA536*k8;
    dwdx[72] = 1.0*APC*k1;
    dwdx[73] = -1.0*VaA3*k8;
    dwdx[74] = -1.0*VaA53*k8;
    dwdx[75] = -1.0*VaA36*k8;
    dwdx[76] = -1.0*VaA536*k7;
    dwdx[77] = -1.0*k2;
    dwdx[78] = -1.0*k2;
    dwdx[79] = 1.0*k3;
    dwdx[80] = 1.0*k5;
}