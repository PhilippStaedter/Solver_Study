#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_brands1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[12] = 1;
    dydx[24] = 1;
    dydx[36] = 1;
    dydx[48] = 1;
    dydx[60] = 1;
    dydx[72] = 1;
    dydx[84] = 1;
    dydx[96] = 1;
    dydx[108] = 1;
    dydx[120] = 1;
}