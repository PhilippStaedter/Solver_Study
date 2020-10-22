#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_nyman3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = ins*k1a + k1aBasic;
    dwdx[1] = k21*(IRiP*k22 + IRp)/(Xp*km23 + 1);
    dwdx[2] = X*k3;
    dwdx[3] = km2;
    dwdx[4] = k1r;
    dwdx[5] = Xp*k1f/(Xp + 1) + k1e;
    dwdx[6] = IRS*k21*k22/(Xp*km23 + 1);
    dwdx[7] = k1b;
    dwdx[8] = k1c;
    dwdx[9] = k1d;
    dwdx[10] = k1g;
    dwdx[11] = IRS*k21/(Xp*km23 + 1);
    dwdx[12] = IRSiP*k3;
    dwdx[13] = IRiP*(-Xp*k1f/pow(Xp + 1, 2) + k1f/(Xp + 1));
    dwdx[14] = -IRS*k21*km23*(IRiP*k22 + IRp)/pow(Xp*km23 + 1, 2);
    dwdx[15] = km3;
}