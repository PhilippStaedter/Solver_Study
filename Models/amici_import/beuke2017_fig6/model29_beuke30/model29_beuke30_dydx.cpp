#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model29_beuke30(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[34] = 1;
    dydx[68] = 1;
    dydx[102] = 1;
    dydx[136] = 1;
    dydx[170] = 1;
    dydx[204] = 1;
    dydx[238] = 1;
    dydx[272] = 1;
    dydx[306] = 1;
    dydx[340] = 1;
    dydx[374] = 1;
    dydx[408] = 1;
    dydx[442] = 1;
    dydx[476] = 1;
    dydx[510] = 1;
    dydx[544] = 1;
    dydx[578] = 1;
    dydx[612] = 1;
    dydx[646] = 1;
    dydx[680] = 1;
    dydx[714] = 1;
    dydx[748] = 1;
    dydx[782] = 1;
    dydx[816] = 1;
    dydx[850] = 1;
    dydx[884] = 1;
    dydx[918] = 1;
    dydx[952] = 1;
    dydx[986] = 1;
    dydx[1020] = 1;
    dydx[1054] = 1;
    dydx[1088] = 1;
}