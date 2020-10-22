#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_panteleev1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[10] = 1;
    dydx[20] = 1;
    dydx[30] = 1;
    dydx[40] = 1;
    dydx[50] = 1;
    dydx[60] = 1;
    dydx[70] = 1;
    dydx[80] = 1;
}