#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Eungdamrong2007(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[47] = 1;
    dydx[94] = 1;
    dydx[141] = 1;
    dydx[188] = 1;
    dydx[235] = 1;
    dydx[282] = 1;
    dydx[329] = 1;
    dydx[376] = 1;
    dydx[423] = 1;
    dydx[470] = 1;
    dydx[517] = 1;
    dydx[564] = 1;
    dydx[611] = 1;
    dydx[658] = 1;
    dydx[705] = 1;
    dydx[752] = 1;
    dydx[799] = 1;
    dydx[846] = 1;
    dydx[893] = 1;
    dydx[940] = 1;
    dydx[987] = 1;
    dydx[1034] = 1;
    dydx[1081] = 1;
    dydx[1128] = 1;
    dydx[1175] = 1;
    dydx[1222] = 1;
    dydx[1269] = 1;
    dydx[1316] = 1;
    dydx[1363] = 1;
    dydx[1410] = 1;
    dydx[1457] = 1;
    dydx[1504] = 1;
    dydx[1551] = 1;
    dydx[1598] = 1;
    dydx[1645] = 1;
    dydx[1692] = 1;
    dydx[1739] = 1;
    dydx[1786] = 1;
    dydx[1833] = 1;
    dydx[1880] = 1;
    dydx[1927] = 1;
    dydx[1974] = 1;
    dydx[2021] = 1;
    dydx[2068] = 1;
    dydx[2115] = 1;
}