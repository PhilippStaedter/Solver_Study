#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Proctor2013a(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[76] = 1;
    dydx[152] = 1;
    dydx[228] = 1;
    dydx[304] = 1;
    dydx[380] = 1;
    dydx[456] = 1;
    dydx[532] = 1;
    dydx[608] = 1;
    dydx[684] = 1;
    dydx[760] = 1;
    dydx[836] = 1;
    dydx[912] = 1;
    dydx[988] = 1;
    dydx[1064] = 1;
    dydx[1140] = 1;
    dydx[1216] = 1;
    dydx[1292] = 1;
    dydx[1368] = 1;
    dydx[1444] = 1;
    dydx[1520] = 1;
    dydx[1596] = 1;
    dydx[1672] = 1;
    dydx[1748] = 1;
    dydx[1824] = 1;
    dydx[1900] = 1;
    dydx[1976] = 1;
    dydx[2052] = 1;
    dydx[2128] = 1;
    dydx[2204] = 1;
    dydx[2280] = 1;
    dydx[2356] = 1;
    dydx[2432] = 1;
    dydx[2508] = 1;
    dydx[2584] = 1;
    dydx[2660] = 1;
    dydx[2736] = 1;
    dydx[2812] = 1;
    dydx[2888] = 1;
    dydx[2964] = 1;
    dydx[3040] = 1;
    dydx[3116] = 1;
    dydx[3192] = 1;
    dydx[3268] = 1;
    dydx[3344] = 1;
    dydx[3420] = 1;
    dydx[3496] = 1;
    dydx[3572] = 1;
    dydx[3648] = 1;
    dydx[3724] = 1;
    dydx[3800] = 1;
    dydx[3876] = 1;
    dydx[3952] = 1;
    dydx[4028] = 1;
    dydx[4104] = 1;
    dydx[4180] = 1;
    dydx[4256] = 1;
    dydx[4332] = 1;
    dydx[4408] = 1;
    dydx[4484] = 1;
    dydx[4560] = 1;
    dydx[4636] = 1;
    dydx[4712] = 1;
    dydx[4788] = 1;
    dydx[4864] = 1;
    dydx[4940] = 1;
    dydx[5016] = 1;
    dydx[5092] = 1;
    dydx[5168] = 1;
    dydx[5244] = 1;
    dydx[5320] = 1;
    dydx[5396] = 1;
    dydx[5472] = 1;
    dydx[5548] = 1;
    dydx[5624] = 1;
}