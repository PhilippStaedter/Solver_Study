#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_aguda1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[13] = 1;
    dydx[26] = 1;
    dydx[39] = 1;
    dydx[52] = 1;
    dydx[65] = 1;
    dydx[78] = 1;
    dydx[91] = 1;
    dydx[104] = 1;
    dydx[117] = 1;
    dydx[130] = 1;
    dydx[143] = 1;
}