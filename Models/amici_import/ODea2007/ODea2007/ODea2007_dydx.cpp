#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_ODea2007(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[25] = 1;
    dydx[50] = 1;
    dydx[75] = 1;
    dydx[100] = 1;
    dydx[125] = 1;
    dydx[150] = 1;
    dydx[175] = 1;
    dydx[200] = 1;
    dydx[225] = 1;
    dydx[250] = 1;
    dydx[275] = 1;
    dydx[300] = 1;
    dydx[325] = 1;
    dydx[350] = 1;
    dydx[375] = 1;
    dydx[400] = 1;
    dydx[425] = 1;
    dydx[450] = 1;
    dydx[475] = 1;
    dydx[500] = 1;
    dydx[525] = 1;
    dydx[550] = 1;
    dydx[575] = 1;
}