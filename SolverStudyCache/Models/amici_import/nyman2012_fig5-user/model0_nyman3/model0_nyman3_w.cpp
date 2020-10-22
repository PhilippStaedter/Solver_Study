#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_nyman3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = IR*ins*k1a + IR*k1aBasic;
    w[1] = IRins*k1b;
    w[2] = IRins*k1c;
    w[3] = IRp*k1d;
    w[4] = IRiP*(Xp*k1f/(Xp + 1) + k1e);
    w[5] = IRp*k1g;
    w[6] = IRi*k1r;
    w[7] = IRS*k21*(IRiP*k22 + IRp)/(Xp*km23 + 1);
    w[8] = IRSiP*X*k3;
    w[9] = IRSiP*km2;
    w[10] = Xp*km3;
}