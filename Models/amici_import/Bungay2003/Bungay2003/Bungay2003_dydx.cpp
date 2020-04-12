#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Bungay2003(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[75] = 1;
    dydx[150] = 1;
    dydx[225] = 1;
    dydx[300] = 1;
    dydx[375] = 1;
    dydx[450] = 1;
    dydx[525] = 1;
    dydx[600] = 1;
    dydx[675] = 1;
    dydx[750] = 1;
    dydx[825] = 1;
    dydx[900] = 1;
    dydx[975] = 1;
    dydx[1050] = 1;
    dydx[1125] = 1;
    dydx[1200] = 1;
    dydx[1275] = 1;
    dydx[1350] = 1;
    dydx[1425] = 1;
    dydx[1500] = 1;
    dydx[1575] = 1;
    dydx[1650] = 1;
    dydx[1725] = 1;
    dydx[1800] = 1;
    dydx[1875] = 1;
    dydx[1950] = 1;
    dydx[2025] = 1;
    dydx[2100] = 1;
    dydx[2175] = 1;
    dydx[2250] = 1;
    dydx[2325] = 1;
    dydx[2400] = 1;
    dydx[2475] = 1;
    dydx[2550] = 1;
    dydx[2625] = 1;
    dydx[2700] = 1;
    dydx[2775] = 1;
    dydx[2850] = 1;
    dydx[2925] = 1;
    dydx[3000] = 1;
    dydx[3075] = 1;
    dydx[3150] = 1;
    dydx[3225] = 1;
    dydx[3300] = 1;
    dydx[3375] = 1;
    dydx[3450] = 1;
    dydx[3525] = 1;
    dydx[3600] = 1;
    dydx[3675] = 1;
    dydx[3750] = 1;
    dydx[3825] = 1;
    dydx[3900] = 1;
    dydx[3975] = 1;
    dydx[4050] = 1;
    dydx[4125] = 1;
    dydx[4200] = 1;
    dydx[4275] = 1;
    dydx[4350] = 1;
    dydx[4425] = 1;
    dydx[4500] = 1;
    dydx[4575] = 1;
    dydx[4650] = 1;
    dydx[4725] = 1;
    dydx[4800] = 1;
    dydx[4875] = 1;
    dydx[4950] = 1;
    dydx[5025] = 1;
    dydx[5100] = 1;
    dydx[5175] = 1;
    dydx[5250] = 1;
    dydx[5325] = 1;
    dydx[5400] = 1;
    dydx[5475] = 1;
}