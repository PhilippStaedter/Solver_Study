#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_kholodenko1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[24] = 1;
    dydx[48] = 1;
    dydx[72] = 1;
    dydx[96] = 1;
    dydx[120] = 1;
    dydx[144] = 1;
    dydx[168] = 1;
    dydx[192] = 1;
    dydx[216] = 1;
    dydx[240] = 1;
    dydx[264] = 1;
    dydx[288] = 1;
    dydx[312] = 1;
    dydx[336] = 1;
    dydx[360] = 1;
    dydx[384] = 1;
    dydx[408] = 1;
    dydx[432] = 1;
    dydx[456] = 1;
    dydx[480] = 1;
    dydx[504] = 1;
    dydx[528] = 1;
}