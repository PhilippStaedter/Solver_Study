#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Bungay2006(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[79] = 1;
    dydx[158] = 1;
    dydx[237] = 1;
    dydx[316] = 1;
    dydx[395] = 1;
    dydx[474] = 1;
    dydx[553] = 1;
    dydx[632] = 1;
    dydx[711] = 1;
    dydx[790] = 1;
    dydx[869] = 1;
    dydx[948] = 1;
    dydx[1027] = 1;
    dydx[1106] = 1;
    dydx[1185] = 1;
    dydx[1264] = 1;
    dydx[1343] = 1;
    dydx[1422] = 1;
    dydx[1501] = 1;
    dydx[1580] = 1;
    dydx[1659] = 1;
    dydx[1738] = 1;
    dydx[1817] = 1;
    dydx[1896] = 1;
    dydx[1975] = 1;
    dydx[2054] = 1;
    dydx[2133] = 1;
    dydx[2212] = 1;
    dydx[2291] = 1;
    dydx[2370] = 1;
    dydx[2449] = 1;
    dydx[2528] = 1;
    dydx[2607] = 1;
    dydx[2686] = 1;
    dydx[2765] = 1;
    dydx[2844] = 1;
    dydx[2923] = 1;
    dydx[3002] = 1;
    dydx[3081] = 1;
    dydx[3160] = 1;
    dydx[3239] = 1;
    dydx[3318] = 1;
    dydx[3397] = 1;
    dydx[3476] = 1;
    dydx[3555] = 1;
    dydx[3634] = 1;
    dydx[3713] = 1;
    dydx[3792] = 1;
    dydx[3871] = 1;
    dydx[3950] = 1;
    dydx[4029] = 1;
    dydx[4108] = 1;
    dydx[4187] = 1;
    dydx[4266] = 1;
    dydx[4345] = 1;
    dydx[4424] = 1;
    dydx[4503] = 1;
    dydx[4582] = 1;
    dydx[4661] = 1;
    dydx[4740] = 1;
    dydx[4819] = 1;
    dydx[4898] = 1;
    dydx[4977] = 1;
    dydx[5056] = 1;
    dydx[5135] = 1;
    dydx[5214] = 1;
    dydx[5293] = 1;
    dydx[5372] = 1;
    dydx[5451] = 1;
    dydx[5530] = 1;
    dydx[5609] = 1;
    dydx[5688] = 1;
    dydx[5767] = 1;
    dydx[5846] = 1;
    dydx[5925] = 1;
    dydx[6004] = 1;
    dydx[6083] = 1;
}