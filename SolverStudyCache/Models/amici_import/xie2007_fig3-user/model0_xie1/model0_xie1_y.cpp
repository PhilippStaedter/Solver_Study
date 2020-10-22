#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_xie1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CC;
    y[1] = CCPT;
    y[2] = CLK;
    y[3] = CYC;
    y[4] = PDP;
    y[5] = PER;
    y[6] = PT;
    y[7] = TIM;
    y[8] = VRI;
    y[9] = clkm;
    y[10] = clkp;
    y[11] = pdpm;
    y[12] = pdpp;
    y[13] = perm;
    y[14] = perp;
    y[15] = prcpdp;
    y[16] = prcper;
    y[17] = prct;
    y[18] = prcv;
    y[19] = prpc;
    y[20] = prvc;
    y[21] = timm;
    y[22] = timp;
    y[23] = vrim;
    y[24] = vrip;
}