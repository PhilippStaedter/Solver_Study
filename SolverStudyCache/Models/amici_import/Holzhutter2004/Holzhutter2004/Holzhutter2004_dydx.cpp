#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Holzhutter2004(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[46] = 1;
    dydx[92] = 1;
    dydx[138] = 1;
    dydx[184] = 1;
    dydx[230] = 1;
    dydx[276] = 1;
    dydx[322] = 1;
    dydx[368] = 1;
    dydx[414] = 1;
    dydx[460] = 1;
    dydx[506] = 1;
    dydx[552] = 1;
    dydx[598] = 1;
    dydx[644] = 1;
    dydx[690] = 1;
    dydx[736] = 1;
    dydx[782] = 1;
    dydx[828] = 1;
    dydx[874] = 1;
    dydx[920] = 1;
    dydx[966] = 1;
    dydx[1012] = 1;
    dydx[1058] = 1;
    dydx[1104] = 1;
    dydx[1150] = 1;
    dydx[1196] = 1;
    dydx[1242] = 1;
    dydx[1288] = 1;
    dydx[1334] = 1;
    dydx[1380] = 1;
    dydx[1426] = 1;
    dydx[1472] = 1;
    dydx[1518] = 1;
    dydx[1564] = 1;
    dydx[1610] = 1;
    dydx[1656] = 1;
    dydx[1702] = 1;
    dydx[1748] = 1;
    dydx[1794] = 1;
    dydx[1840] = 1;
    dydx[1886] = 1;
    dydx[1932] = 1;
    dydx[1978] = 1;
    dydx[2024] = 1;
}