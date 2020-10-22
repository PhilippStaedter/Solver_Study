#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Levchenko2000a(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[87] = 1;
    dydx[174] = 1;
    dydx[261] = 1;
    dydx[348] = 1;
    dydx[435] = 1;
    dydx[522] = 1;
    dydx[609] = 1;
    dydx[696] = 1;
    dydx[783] = 1;
    dydx[870] = 1;
    dydx[957] = 1;
    dydx[1044] = 1;
    dydx[1131] = 1;
    dydx[1218] = 1;
    dydx[1305] = 1;
    dydx[1392] = 1;
    dydx[1479] = 1;
    dydx[1566] = 1;
    dydx[1653] = 1;
    dydx[1740] = 1;
    dydx[1827] = 1;
    dydx[1914] = 1;
    dydx[2001] = 1;
    dydx[2088] = 1;
    dydx[2175] = 1;
    dydx[2262] = 1;
    dydx[2349] = 1;
    dydx[2436] = 1;
    dydx[2523] = 1;
    dydx[2610] = 1;
    dydx[2697] = 1;
    dydx[2784] = 1;
    dydx[2871] = 1;
    dydx[2958] = 1;
    dydx[3045] = 1;
    dydx[3132] = 1;
    dydx[3219] = 1;
    dydx[3306] = 1;
    dydx[3393] = 1;
    dydx[3480] = 1;
    dydx[3567] = 1;
    dydx[3654] = 1;
    dydx[3741] = 1;
    dydx[3828] = 1;
    dydx[3915] = 1;
    dydx[4002] = 1;
    dydx[4089] = 1;
    dydx[4176] = 1;
    dydx[4263] = 1;
    dydx[4350] = 1;
    dydx[4437] = 1;
    dydx[4524] = 1;
    dydx[4611] = 1;
    dydx[4698] = 1;
    dydx[4785] = 1;
    dydx[4872] = 1;
    dydx[4959] = 1;
    dydx[5046] = 1;
    dydx[5133] = 1;
    dydx[5220] = 1;
    dydx[5307] = 1;
    dydx[5394] = 1;
    dydx[5481] = 1;
    dydx[5568] = 1;
    dydx[5655] = 1;
    dydx[5742] = 1;
    dydx[5829] = 1;
    dydx[5916] = 1;
    dydx[6003] = 1;
    dydx[6090] = 1;
    dydx[6177] = 1;
    dydx[6264] = 1;
    dydx[6351] = 1;
    dydx[6438] = 1;
    dydx[6525] = 1;
    dydx[6612] = 1;
    dydx[6699] = 1;
    dydx[6786] = 1;
    dydx[6873] = 1;
    dydx[6960] = 1;
    dydx[7047] = 1;
    dydx[7134] = 1;
    dydx[7221] = 1;
    dydx[7308] = 1;
    dydx[7395] = 1;
}