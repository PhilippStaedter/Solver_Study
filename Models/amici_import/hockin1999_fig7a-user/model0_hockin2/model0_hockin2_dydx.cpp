#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_hockin2(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[31] = 1;
    dydx[62] = 1;
    dydx[93] = 1;
    dydx[124] = 1;
    dydx[155] = 1;
    dydx[186] = 1;
    dydx[217] = 1;
    dydx[248] = 1;
    dydx[279] = 1;
    dydx[310] = 1;
    dydx[341] = 1;
    dydx[372] = 1;
    dydx[403] = 1;
    dydx[434] = 1;
    dydx[465] = 1;
    dydx[496] = 1;
    dydx[527] = 1;
    dydx[558] = 1;
    dydx[589] = 1;
    dydx[620] = 1;
    dydx[651] = 1;
    dydx[682] = 1;
    dydx[713] = 1;
    dydx[744] = 1;
    dydx[775] = 1;
    dydx[806] = 1;
    dydx[837] = 1;
    dydx[868] = 1;
    dydx[899] = 1;
}