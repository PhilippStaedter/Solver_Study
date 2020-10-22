#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_kofahl(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[38] = 1;
    dydx[76] = 1;
    dydx[114] = 1;
    dydx[152] = 1;
    dydx[190] = 1;
    dydx[228] = 1;
    dydx[266] = 1;
    dydx[304] = 1;
    dydx[342] = 1;
    dydx[380] = 1;
    dydx[418] = 1;
    dydx[456] = 1;
    dydx[494] = 1;
    dydx[532] = 1;
    dydx[570] = 1;
    dydx[608] = 1;
    dydx[646] = 1;
    dydx[684] = 1;
    dydx[722] = 1;
    dydx[760] = 1;
    dydx[798] = 1;
    dydx[836] = 1;
    dydx[874] = 1;
    dydx[912] = 1;
    dydx[950] = 1;
    dydx[988] = 1;
    dydx[1026] = 1;
    dydx[1064] = 1;
    dydx[1102] = 1;
    dydx[1140] = 1;
    dydx[1178] = 1;
    dydx[1216] = 1;
    dydx[1254] = 1;
    dydx[1292] = 1;
    dydx[1330] = 1;
    dydx[1368] = 1;
}