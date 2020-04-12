#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model6_levchenko1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[23] = 1;
    dydx[46] = 1;
    dydx[69] = 1;
    dydx[92] = 1;
    dydx[115] = 1;
    dydx[138] = 1;
    dydx[161] = 1;
    dydx[184] = 1;
    dydx[207] = 1;
    dydx[230] = 1;
    dydx[253] = 1;
    dydx[276] = 1;
    dydx[299] = 1;
    dydx[322] = 1;
    dydx[345] = 1;
    dydx[368] = 1;
    dydx[391] = 1;
    dydx[414] = 1;
    dydx[437] = 1;
    dydx[460] = 1;
    dydx[483] = 1;
}