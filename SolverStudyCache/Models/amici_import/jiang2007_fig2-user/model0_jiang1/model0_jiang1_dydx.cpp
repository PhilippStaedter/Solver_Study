#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_jiang1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[60] = 1;
    dydx[120] = 1;
    dydx[180] = 1;
    dydx[240] = 1;
    dydx[300] = 1;
    dydx[360] = 1;
    dydx[420] = 1;
    dydx[480] = 1;
    dydx[540] = 1;
    dydx[600] = 1;
    dydx[660] = 1;
    dydx[720] = 1;
    dydx[780] = 1;
    dydx[840] = 1;
    dydx[900] = 1;
    dydx[960] = 1;
    dydx[1020] = 1;
    dydx[1080] = 1;
    dydx[1140] = 1;
    dydx[1200] = 1;
    dydx[1260] = 1;
    dydx[1320] = 1;
    dydx[1380] = 1;
    dydx[1440] = 1;
    dydx[1500] = 1;
    dydx[1560] = 1;
    dydx[1620] = 1;
    dydx[1680] = 1;
    dydx[1740] = 1;
    dydx[1800] = 1;
    dydx[1860] = 1;
    dydx[1920] = 1;
    dydx[1980] = 1;
    dydx[2040] = 1;
    dydx[2100] = 1;
    dydx[2160] = 1;
    dydx[2220] = 1;
    dydx[2280] = 1;
    dydx[2340] = 1;
    dydx[2400] = 1;
    dydx[2460] = 1;
    dydx[2520] = 1;
    dydx[2580] = 1;
    dydx[2640] = 1;
    dydx[2700] = 1;
    dydx[2760] = 1;
    dydx[2820] = 1;
    dydx[2880] = 1;
    dydx[2940] = 1;
    dydx[3000] = 1;
    dydx[3060] = 1;
    dydx[3120] = 1;
    dydx[3180] = 1;
    dydx[3240] = 1;
    dydx[3300] = 1;
    dydx[3360] = 1;
    dydx[3420] = 1;
    dydx[3480] = 1;
}