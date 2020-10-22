#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model5_mcauley1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[35] = 1;
    dydx[70] = 1;
    dydx[105] = 1;
    dydx[140] = 1;
    dydx[175] = 1;
    dydx[210] = 1;
    dydx[245] = 1;
    dydx[280] = 1;
    dydx[315] = 1;
    dydx[350] = 1;
    dydx[385] = 1;
    dydx[420] = 1;
    dydx[455] = 1;
    dydx[490] = 1;
    dydx[525] = 1;
    dydx[560] = 1;
    dydx[595] = 1;
    dydx[630] = 1;
    dydx[665] = 1;
    dydx[700] = 1;
    dydx[735] = 1;
    dydx[770] = 1;
    dydx[805] = 1;
    dydx[840] = 1;
    dydx[875] = 1;
    dydx[910] = 1;
    dydx[945] = 1;
    dydx[980] = 1;
    dydx[1015] = 1;
    dydx[1050] = 1;
    dydx[1085] = 1;
    dydx[1120] = 1;
    dydx[1155] = 1;
}