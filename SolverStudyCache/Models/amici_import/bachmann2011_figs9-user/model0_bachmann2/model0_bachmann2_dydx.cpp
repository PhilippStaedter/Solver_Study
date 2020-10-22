#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_bachmann2(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[27] = 1;
    dydx[54] = 1;
    dydx[81] = 1;
    dydx[108] = 1;
    dydx[135] = 1;
    dydx[162] = 1;
    dydx[189] = 1;
    dydx[216] = 1;
    dydx[243] = 1;
    dydx[270] = 1;
    dydx[297] = 1;
    dydx[324] = 1;
    dydx[351] = 1;
    dydx[378] = 1;
    dydx[405] = 1;
    dydx[432] = 1;
    dydx[459] = 1;
    dydx[486] = 1;
    dydx[513] = 1;
    dydx[540] = 1;
    dydx[567] = 1;
    dydx[594] = 1;
    dydx[621] = 1;
    dydx[648] = 1;
    dydx[675] = 1;
}