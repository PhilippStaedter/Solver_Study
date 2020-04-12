#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_chassagnole2(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[19] = 1;
    dydx[38] = 1;
    dydx[57] = 1;
    dydx[76] = 1;
    dydx[95] = 1;
    dydx[114] = 1;
    dydx[133] = 1;
    dydx[152] = 1;
    dydx[171] = 1;
    dydx[190] = 1;
    dydx[209] = 1;
    dydx[228] = 1;
    dydx[247] = 1;
    dydx[266] = 1;
    dydx[285] = 1;
    dydx[304] = 1;
    dydx[323] = 1;
}