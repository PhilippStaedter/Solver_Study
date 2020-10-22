#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_hockin2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*APC*Va*k1 - 1.0*Va_APC*k2;
    w[1] = 1.0*Va56_APC*k3;
    w[2] = 1.0*Va5_APC*k3;
    w[3] = 1.0*Va53_APC*k6;
    w[4] = 1.0*Va_APC*k3;
    w[5] = 1.0*Va3_APC*k6;
    w[6] = 1.0*Va36_APC*k5;
    w[7] = 1.0*Va3_APC*k5;
    w[8] = -1.0*HC*LC*k10 + 1.0*Va*k9;
    w[9] = -1.0*HC5*LC*k10 + 1.0*Va5*k9;
    w[10] = -1.0*HC3*LC*k10 + 1.0*Va3*k9;
    w[11] = 1.0*APC*Va3*k1 - 1.0*Va3_APC*k2;
    w[12] = -1.0*HC56*LC*k10 + 1.0*Va56*k9;
    w[13] = -1.0*HC53*LC*k10 + 1.0*Va53*k9;
    w[14] = -1.0*HC36*LC*k10 + 1.0*Va36*k9;
    w[15] = -1.0*HC536*LC*k10 + 1.0*Va536*k9;
    w[16] = 1.0*APC*LC*k1 - 1.0*LC_APC*k2;
    w[17] = 1.0*Va3*k7 - 1.0*VaA3*VaLCA1*k8;
    w[18] = 1.0*Va53*k7 - 1.0*VaA53*VaLCA1*k8;
    w[19] = 1.0*Va36*k7 - 1.0*VaA36*VaLCA1*k8;
    w[20] = 1.0*Va536*k7 - 1.0*VaA536*VaLCA1*k8;
    w[21] = 1.0*Va3_APC*k7 - 1.0*VaA3*VaLCA1_APC*k8;
    w[22] = 1.0*APC*Va5*k1 - 1.0*Va5_APC*k2;
    w[23] = 1.0*Va53_APC*k7 - 1.0*VaA53*VaLCA1_APC*k8;
    w[24] = 1.0*Va36_APC*k7 - 1.0*VaA36*VaLCA1_APC*k8;
    w[25] = 1.0*Va536_APC*k7 - 1.0*VaA536*VaLCA1_APC*k7;
    w[26] = 1.0*APC*VaLCA1*k1 - 1.0*VaLCA1_APC*k2;
    w[27] = 1.0*APC*Va53*k1 - 1.0*Va53_APC*k2;
    w[28] = 1.0*APC*Va56*k1 - 1.0*Va56_APC*k2;
    w[29] = 1.0*APC*Va36*k1 - 1.0*Va36_APC*k2;
    w[30] = 1.0*APC*Va536*k1 - 1.0*Va536_APC*k2;
    w[31] = 1.0*Va_APC*k5;
    w[32] = 1.0*Va5_APC*k6;
}