#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_lee4(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[15] = 1;
    dydx[30] = 1;
    dydx[45] = 1;
    dydx[60] = 1;
    dydx[75] = 1;
    dydx[90] = 1;
    dydx[105] = 1;
    dydx[120] = 1;
    dydx[135] = 1;
    dydx[150] = 1;
    dydx[165] = 1;
    dydx[180] = 1;
    dydx[195] = 1;
}