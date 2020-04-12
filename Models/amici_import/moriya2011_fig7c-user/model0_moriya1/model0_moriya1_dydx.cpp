#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_model0_moriya1(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[50] = 1;
    dydx[100] = 1;
    dydx[150] = 1;
    dydx[200] = 1;
    dydx[250] = 1;
    dydx[300] = 1;
    dydx[350] = 1;
    dydx[400] = 1;
    dydx[450] = 1;
    dydx[500] = 1;
    dydx[550] = 1;
    dydx[600] = 1;
    dydx[650] = 1;
    dydx[700] = 1;
    dydx[750] = 1;
    dydx[800] = 1;
    dydx[850] = 1;
    dydx[900] = 1;
    dydx[950] = 1;
    dydx[1000] = 1;
    dydx[1050] = 1;
    dydx[1100] = 1;
    dydx[1150] = 1;
    dydx[1200] = 1;
    dydx[1250] = 1;
    dydx[1300] = 1;
    dydx[1350] = 1;
    dydx[1400] = 1;
    dydx[1450] = 1;
    dydx[1500] = 1;
    dydx[1550] = 1;
    dydx[1600] = 1;
    dydx[1650] = 1;
    dydx[1700] = 1;
    dydx[1750] = 1;
    dydx[1800] = 1;
    dydx[1850] = 1;
    dydx[1900] = 1;
    dydx[1950] = 1;
    dydx[2000] = 1;
    dydx[2050] = 1;
    dydx[2100] = 1;
    dydx[2150] = 1;
    dydx[2200] = 1;
    dydx[2250] = 1;
    dydx[2300] = 1;
    dydx[2350] = 1;
    dydx[2400] = 1;
}