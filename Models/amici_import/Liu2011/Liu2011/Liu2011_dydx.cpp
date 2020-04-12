#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Liu2011(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[43] = 1;
    dydx[86] = 1;
    dydx[129] = 1;
    dydx[172] = 1;
    dydx[215] = 1;
    dydx[258] = 1;
    dydx[301] = 1;
    dydx[344] = 1;
    dydx[387] = 1;
    dydx[430] = 1;
    dydx[473] = 1;
    dydx[516] = 1;
    dydx[559] = 1;
    dydx[602] = 1;
    dydx[645] = 1;
    dydx[688] = 1;
    dydx[731] = 1;
    dydx[774] = 1;
    dydx[817] = 1;
    dydx[860] = 1;
    dydx[903] = 1;
    dydx[946] = 1;
    dydx[989] = 1;
    dydx[1032] = 1;
    dydx[1075] = 1;
    dydx[1118] = 1;
    dydx[1161] = 1;
    dydx[1204] = 1;
    dydx[1247] = 1;
    dydx[1290] = 1;
    dydx[1333] = 1;
    dydx[1376] = 1;
    dydx[1419] = 1;
    dydx[1462] = 1;
    dydx[1505] = 1;
    dydx[1548] = 1;
    dydx[1591] = 1;
    dydx[1634] = 1;
    dydx[1677] = 1;
    dydx[1720] = 1;
    dydx[1763] = 1;
}