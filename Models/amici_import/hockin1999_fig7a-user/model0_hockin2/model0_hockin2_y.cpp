#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_hockin2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = APC;
    y[1] = HC;
    y[2] = HC3;
    y[3] = HC36;
    y[4] = HC5;
    y[5] = HC53;
    y[6] = HC536;
    y[7] = HC56;
    y[8] = LC;
    y[9] = LC_APC;
    y[10] = Va;
    y[11] = Va3;
    y[12] = Va36;
    y[13] = Va36_APC;
    y[14] = Va3_APC;
    y[15] = Va5;
    y[16] = Va53;
    y[17] = Va536;
    y[18] = Va536_APC;
    y[19] = Va53_APC;
    y[20] = Va56;
    y[21] = Va56_APC;
    y[22] = Va5_APC;
    y[23] = VaA3;
    y[24] = VaA36;
    y[25] = VaA53;
    y[26] = VaA536;
    y[27] = VaLCA1;
    y[28] = VaLCA1_APC;
    y[29] = Va_APC;
}