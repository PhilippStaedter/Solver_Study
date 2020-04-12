#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_zi1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[17] = 1;
    dydx[34] = 1;
    dydx[51] = 1;
    dydx[68] = 1;
    dydx[85] = 1;
    dydx[102] = 1;
    dydx[119] = 1;
    dydx[136] = 1;
    dydx[153] = 1;
    dydx[170] = 1;
    dydx[187] = 1;
    dydx[204] = 1;
    dydx[221] = 1;
    dydx[238] = 1;
    dydx[255] = 1;
}