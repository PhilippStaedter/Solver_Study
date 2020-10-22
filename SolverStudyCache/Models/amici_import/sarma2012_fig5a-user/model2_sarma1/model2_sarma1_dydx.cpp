#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model2_sarma1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[28] = 1;
    dydx[56] = 1;
    dydx[84] = 1;
    dydx[112] = 1;
    dydx[140] = 1;
    dydx[168] = 1;
    dydx[196] = 1;
    dydx[224] = 1;
    dydx[252] = 1;
    dydx[280] = 1;
    dydx[308] = 1;
    dydx[336] = 1;
    dydx[364] = 1;
    dydx[392] = 1;
    dydx[420] = 1;
    dydx[448] = 1;
    dydx[476] = 1;
    dydx[504] = 1;
    dydx[532] = 1;
    dydx[560] = 1;
    dydx[588] = 1;
    dydx[616] = 1;
    dydx[644] = 1;
    dydx[672] = 1;
    dydx[700] = 1;
    dydx[728] = 1;
}