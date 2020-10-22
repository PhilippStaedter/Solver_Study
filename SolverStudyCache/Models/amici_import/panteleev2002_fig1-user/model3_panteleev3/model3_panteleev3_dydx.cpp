#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model3_panteleev3(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[9] = 1;
    dydx[18] = 1;
    dydx[27] = 1;
    dydx[36] = 1;
    dydx[45] = 1;
    dydx[54] = 1;
    dydx[63] = 1;
}