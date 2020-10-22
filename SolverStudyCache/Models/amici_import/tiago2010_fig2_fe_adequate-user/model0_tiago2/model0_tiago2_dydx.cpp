#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_tiago2(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[18] = 1;
    dydx[36] = 1;
    dydx[54] = 1;
    dydx[72] = 1;
    dydx[90] = 1;
    dydx[108] = 1;
    dydx[126] = 1;
    dydx[144] = 1;
    dydx[162] = 1;
    dydx[180] = 1;
    dydx[198] = 1;
    dydx[216] = 1;
    dydx[234] = 1;
    dydx[252] = 1;
    dydx[270] = 1;
    dydx[288] = 1;
}