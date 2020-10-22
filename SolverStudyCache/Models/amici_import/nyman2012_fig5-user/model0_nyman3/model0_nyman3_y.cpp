#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_nyman3(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IR;
    y[1] = IRS;
    y[2] = IRSiP;
    y[3] = IRi;
    y[4] = IRiP;
    y[5] = IRins;
    y[6] = IRp;
    y[7] = X;
    y[8] = Xp;
}