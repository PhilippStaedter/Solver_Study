#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model3_sarma7(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[14] = 1;
    dydx[28] = 1;
    dydx[42] = 1;
    dydx[56] = 1;
    dydx[70] = 1;
    dydx[84] = 1;
    dydx[98] = 1;
    dydx[112] = 1;
    dydx[126] = 1;
    dydx[140] = 1;
    dydx[154] = 1;
    dydx[168] = 1;
}