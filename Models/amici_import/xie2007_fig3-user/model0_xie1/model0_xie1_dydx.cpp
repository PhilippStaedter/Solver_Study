#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_xie1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[26] = 1;
    dydx[52] = 1;
    dydx[78] = 1;
    dydx[104] = 1;
    dydx[130] = 1;
    dydx[156] = 1;
    dydx[182] = 1;
    dydx[208] = 1;
    dydx[234] = 1;
    dydx[260] = 1;
    dydx[286] = 1;
    dydx[312] = 1;
    dydx[338] = 1;
    dydx[364] = 1;
    dydx[390] = 1;
    dydx[416] = 1;
    dydx[442] = 1;
    dydx[468] = 1;
    dydx[494] = 1;
    dydx[520] = 1;
    dydx[546] = 1;
    dydx[572] = 1;
    dydx[598] = 1;
    dydx[624] = 1;
}