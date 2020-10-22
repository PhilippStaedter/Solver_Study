#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_levchenko2(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[32] = 1;
    dydx[64] = 1;
    dydx[96] = 1;
    dydx[128] = 1;
    dydx[160] = 1;
    dydx[192] = 1;
    dydx[224] = 1;
    dydx[256] = 1;
    dydx[288] = 1;
    dydx[320] = 1;
    dydx[352] = 1;
    dydx[384] = 1;
    dydx[416] = 1;
    dydx[448] = 1;
    dydx[480] = 1;
    dydx[512] = 1;
    dydx[544] = 1;
    dydx[576] = 1;
    dydx[608] = 1;
    dydx[640] = 1;
    dydx[672] = 1;
    dydx[704] = 1;
    dydx[736] = 1;
    dydx[768] = 1;
    dydx[800] = 1;
    dydx[832] = 1;
    dydx[864] = 1;
    dydx[896] = 1;
    dydx[928] = 1;
    dydx[960] = 1;
}