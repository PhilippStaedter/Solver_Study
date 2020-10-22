#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Bungay2006a(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[55] = 1;
    dydx[110] = 1;
    dydx[165] = 1;
    dydx[220] = 1;
    dydx[275] = 1;
    dydx[330] = 1;
    dydx[385] = 1;
    dydx[440] = 1;
    dydx[495] = 1;
    dydx[550] = 1;
    dydx[605] = 1;
    dydx[660] = 1;
    dydx[715] = 1;
    dydx[770] = 1;
    dydx[825] = 1;
    dydx[880] = 1;
    dydx[935] = 1;
    dydx[990] = 1;
    dydx[1045] = 1;
    dydx[1100] = 1;
    dydx[1155] = 1;
    dydx[1210] = 1;
    dydx[1265] = 1;
    dydx[1320] = 1;
    dydx[1375] = 1;
    dydx[1430] = 1;
    dydx[1485] = 1;
    dydx[1540] = 1;
    dydx[1595] = 1;
    dydx[1650] = 1;
    dydx[1705] = 1;
    dydx[1760] = 1;
    dydx[1815] = 1;
    dydx[1870] = 1;
    dydx[1925] = 1;
    dydx[1980] = 1;
    dydx[2035] = 1;
    dydx[2090] = 1;
    dydx[2145] = 1;
    dydx[2200] = 1;
    dydx[2255] = 1;
    dydx[2310] = 1;
    dydx[2365] = 1;
    dydx[2420] = 1;
    dydx[2475] = 1;
    dydx[2530] = 1;
    dydx[2585] = 1;
    dydx[2640] = 1;
    dydx[2695] = 1;
    dydx[2750] = 1;
    dydx[2805] = 1;
    dydx[2860] = 1;
    dydx[2915] = 1;
}