#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_hatakeyama1_Fig5G(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[37] = 1;
    dydx[74] = 1;
    dydx[111] = 1;
    dydx[148] = 1;
    dydx[185] = 1;
    dydx[222] = 1;
    dydx[259] = 1;
    dydx[296] = 1;
    dydx[333] = 1;
    dydx[370] = 1;
    dydx[407] = 1;
    dydx[444] = 1;
    dydx[481] = 1;
    dydx[518] = 1;
    dydx[555] = 1;
    dydx[592] = 1;
    dydx[629] = 1;
    dydx[666] = 1;
    dydx[703] = 1;
    dydx[740] = 1;
    dydx[777] = 1;
    dydx[814] = 1;
    dydx[851] = 1;
    dydx[888] = 1;
    dydx[925] = 1;
    dydx[962] = 1;
    dydx[999] = 1;
    dydx[1036] = 1;
    dydx[1073] = 1;
    dydx[1110] = 1;
    dydx[1147] = 1;
    dydx[1184] = 1;
    dydx[1221] = 1;
    dydx[1258] = 1;
    dydx[1295] = 1;
}