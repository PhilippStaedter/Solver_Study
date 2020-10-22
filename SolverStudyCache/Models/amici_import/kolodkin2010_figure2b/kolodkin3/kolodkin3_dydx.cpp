#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_kolodkin3(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[11] = 1;
    dydx[22] = 1;
    dydx[33] = 1;
    dydx[44] = 1;
    dydx[55] = 1;
    dydx[66] = 1;
    dydx[77] = 1;
    dydx[88] = 1;
    dydx[99] = 1;
}