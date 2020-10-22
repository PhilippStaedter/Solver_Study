#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Pathak2013(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[58] = 1;
    dydx[116] = 1;
    dydx[174] = 1;
    dydx[232] = 1;
    dydx[290] = 1;
    dydx[348] = 1;
    dydx[406] = 1;
    dydx[464] = 1;
    dydx[522] = 1;
    dydx[580] = 1;
    dydx[638] = 1;
    dydx[696] = 1;
    dydx[754] = 1;
    dydx[812] = 1;
    dydx[870] = 1;
    dydx[928] = 1;
    dydx[986] = 1;
    dydx[1044] = 1;
    dydx[1102] = 1;
    dydx[1160] = 1;
    dydx[1218] = 1;
    dydx[1276] = 1;
    dydx[1334] = 1;
    dydx[1392] = 1;
    dydx[1450] = 1;
    dydx[1508] = 1;
    dydx[1566] = 1;
    dydx[1624] = 1;
    dydx[1682] = 1;
    dydx[1740] = 1;
    dydx[1798] = 1;
    dydx[1856] = 1;
    dydx[1914] = 1;
    dydx[1972] = 1;
    dydx[2030] = 1;
    dydx[2088] = 1;
    dydx[2146] = 1;
    dydx[2204] = 1;
    dydx[2262] = 1;
    dydx[2320] = 1;
    dydx[2378] = 1;
    dydx[2436] = 1;
    dydx[2494] = 1;
    dydx[2552] = 1;
    dydx[2610] = 1;
    dydx[2668] = 1;
    dydx[2726] = 1;
    dydx[2784] = 1;
    dydx[2842] = 1;
    dydx[2900] = 1;
    dydx[2958] = 1;
    dydx[3016] = 1;
    dydx[3074] = 1;
    dydx[3132] = 1;
    dydx[3190] = 1;
    dydx[3248] = 1;
}