#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Singh2006(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[69] = 1;
    dydx[138] = 1;
    dydx[207] = 1;
    dydx[276] = 1;
    dydx[345] = 1;
    dydx[414] = 1;
    dydx[483] = 1;
    dydx[552] = 1;
    dydx[621] = 1;
    dydx[690] = 1;
    dydx[759] = 1;
    dydx[828] = 1;
    dydx[897] = 1;
    dydx[966] = 1;
    dydx[1035] = 1;
    dydx[1104] = 1;
    dydx[1173] = 1;
    dydx[1242] = 1;
    dydx[1311] = 1;
    dydx[1380] = 1;
    dydx[1449] = 1;
    dydx[1518] = 1;
    dydx[1587] = 1;
    dydx[1656] = 1;
    dydx[1725] = 1;
    dydx[1794] = 1;
    dydx[1863] = 1;
    dydx[1932] = 1;
    dydx[2001] = 1;
    dydx[2070] = 1;
    dydx[2139] = 1;
    dydx[2208] = 1;
    dydx[2277] = 1;
    dydx[2346] = 1;
    dydx[2415] = 1;
    dydx[2484] = 1;
    dydx[2553] = 1;
    dydx[2622] = 1;
    dydx[2691] = 1;
    dydx[2760] = 1;
    dydx[2829] = 1;
    dydx[2898] = 1;
    dydx[2967] = 1;
    dydx[3036] = 1;
    dydx[3105] = 1;
    dydx[3174] = 1;
    dydx[3243] = 1;
    dydx[3312] = 1;
    dydx[3381] = 1;
    dydx[3450] = 1;
    dydx[3519] = 1;
    dydx[3588] = 1;
    dydx[3657] = 1;
    dydx[3726] = 1;
    dydx[3795] = 1;
    dydx[3864] = 1;
    dydx[3933] = 1;
    dydx[4002] = 1;
    dydx[4071] = 1;
    dydx[4140] = 1;
    dydx[4209] = 1;
    dydx[4278] = 1;
    dydx[4347] = 1;
    dydx[4416] = 1;
    dydx[4485] = 1;
    dydx[4554] = 1;
    dydx[4623] = 1;
}