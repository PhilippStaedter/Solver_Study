#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Hui2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[66] = 1;
    dydx[132] = 1;
    dydx[198] = 1;
    dydx[264] = 1;
    dydx[330] = 1;
    dydx[396] = 1;
    dydx[462] = 1;
    dydx[528] = 1;
    dydx[594] = 1;
    dydx[660] = 1;
    dydx[726] = 1;
    dydx[792] = 1;
    dydx[858] = 1;
    dydx[924] = 1;
    dydx[990] = 1;
    dydx[1056] = 1;
    dydx[1122] = 1;
    dydx[1188] = 1;
    dydx[1254] = 1;
    dydx[1320] = 1;
    dydx[1386] = 1;
    dydx[1452] = 1;
    dydx[1518] = 1;
    dydx[1584] = 1;
    dydx[1650] = 1;
    dydx[1716] = 1;
    dydx[1782] = 1;
    dydx[1848] = 1;
    dydx[1914] = 1;
    dydx[1980] = 1;
    dydx[2046] = 1;
    dydx[2112] = 1;
    dydx[2178] = 1;
    dydx[2244] = 1;
    dydx[2310] = 1;
    dydx[2376] = 1;
    dydx[2442] = 1;
    dydx[2508] = 1;
    dydx[2574] = 1;
    dydx[2640] = 1;
    dydx[2706] = 1;
    dydx[2772] = 1;
    dydx[2838] = 1;
    dydx[2904] = 1;
    dydx[2970] = 1;
    dydx[3036] = 1;
    dydx[3102] = 1;
    dydx[3168] = 1;
    dydx[3234] = 1;
    dydx[3300] = 1;
    dydx[3366] = 1;
    dydx[3432] = 1;
    dydx[3498] = 1;
    dydx[3564] = 1;
    dydx[3630] = 1;
    dydx[3696] = 1;
    dydx[3762] = 1;
    dydx[3828] = 1;
    dydx[3894] = 1;
    dydx[3960] = 1;
    dydx[4026] = 1;
    dydx[4092] = 1;
    dydx[4158] = 1;
    dydx[4224] = 1;
}