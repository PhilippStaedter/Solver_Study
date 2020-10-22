#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_hockin2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*APC*Va;
            dwdp[11] = 1.0*APC*Va3;
            dwdp[16] = 1.0*APC*LC;
            dwdp[22] = 1.0*APC*Va5;
            dwdp[26] = 1.0*APC*VaLCA1;
            dwdp[27] = 1.0*APC*Va53;
            dwdp[28] = 1.0*APC*Va56;
            dwdp[29] = 1.0*APC*Va36;
            dwdp[30] = 1.0*APC*Va536;
            break;
        case 1:
            dwdp[8] = -1.0*HC*LC;
            dwdp[9] = -1.0*HC5*LC;
            dwdp[10] = -1.0*HC3*LC;
            dwdp[12] = -1.0*HC56*LC;
            dwdp[13] = -1.0*HC53*LC;
            dwdp[14] = -1.0*HC36*LC;
            dwdp[15] = -1.0*HC536*LC;
            break;
        case 2:
            dwdp[0] = -1.0*Va_APC;
            dwdp[11] = -1.0*Va3_APC;
            dwdp[16] = -1.0*LC_APC;
            dwdp[22] = -1.0*Va5_APC;
            dwdp[26] = -1.0*VaLCA1_APC;
            dwdp[27] = -1.0*Va53_APC;
            dwdp[28] = -1.0*Va56_APC;
            dwdp[29] = -1.0*Va36_APC;
            dwdp[30] = -1.0*Va536_APC;
            break;
        case 3:
            dwdp[1] = 1.0*Va56_APC;
            dwdp[2] = 1.0*Va5_APC;
            dwdp[4] = 1.0*Va_APC;
            break;
        case 4:
            dwdp[6] = 1.0*Va36_APC;
            dwdp[7] = 1.0*Va3_APC;
            dwdp[31] = 1.0*Va_APC;
            break;
        case 5:
            dwdp[3] = 1.0*Va53_APC;
            dwdp[5] = 1.0*Va3_APC;
            dwdp[32] = 1.0*Va5_APC;
            break;
        case 6:
            dwdp[17] = 1.0*Va3;
            dwdp[18] = 1.0*Va53;
            dwdp[19] = 1.0*Va36;
            dwdp[20] = 1.0*Va536;
            dwdp[21] = 1.0*Va3_APC;
            dwdp[23] = 1.0*Va53_APC;
            dwdp[24] = 1.0*Va36_APC;
            dwdp[25] = 1.0*Va536_APC - 1.0*VaA536*VaLCA1_APC;
            break;
        case 7:
            dwdp[17] = -1.0*VaA3*VaLCA1;
            dwdp[18] = -1.0*VaA53*VaLCA1;
            dwdp[19] = -1.0*VaA36*VaLCA1;
            dwdp[20] = -1.0*VaA536*VaLCA1;
            dwdp[21] = -1.0*VaA3*VaLCA1_APC;
            dwdp[23] = -1.0*VaA53*VaLCA1_APC;
            dwdp[24] = -1.0*VaA36*VaLCA1_APC;
            break;
        case 8:
            dwdp[8] = 1.0*Va;
            dwdp[9] = 1.0*Va5;
            dwdp[10] = 1.0*Va3;
            dwdp[12] = 1.0*Va56;
            dwdp[13] = 1.0*Va53;
            dwdp[14] = 1.0*Va36;
            dwdp[15] = 1.0*Va536;
            break;
    }
}