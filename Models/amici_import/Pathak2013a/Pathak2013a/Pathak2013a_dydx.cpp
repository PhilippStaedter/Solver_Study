#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Pathak2013a(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[53] = 1;
    dydx[106] = 1;
    dydx[159] = 1;
    dydx[212] = 1;
    dydx[265] = 1;
    dydx[318] = 1;
    dydx[371] = 1;
    dydx[424] = 1;
    dydx[477] = 1;
    dydx[530] = 1;
    dydx[583] = 1;
    dydx[636] = 1;
    dydx[689] = 1;
    dydx[742] = 1;
    dydx[795] = 1;
    dydx[848] = 1;
    dydx[901] = 1;
    dydx[954] = 1;
    dydx[1007] = 1;
    dydx[1060] = 1;
    dydx[1113] = 1;
    dydx[1166] = 1;
    dydx[1219] = 1;
    dydx[1272] = 1;
    dydx[1325] = 1;
    dydx[1378] = 1;
    dydx[1431] = 1;
    dydx[1484] = 1;
    dydx[1537] = 1;
    dydx[1590] = 1;
    dydx[1643] = 1;
    dydx[1696] = 1;
    dydx[1749] = 1;
    dydx[1802] = 1;
    dydx[1855] = 1;
    dydx[1908] = 1;
    dydx[1961] = 1;
    dydx[2014] = 1;
    dydx[2067] = 1;
    dydx[2120] = 1;
    dydx[2173] = 1;
    dydx[2226] = 1;
    dydx[2279] = 1;
    dydx[2332] = 1;
    dydx[2385] = 1;
    dydx[2438] = 1;
    dydx[2491] = 1;
    dydx[2544] = 1;
    dydx[2597] = 1;
    dydx[2650] = 1;
    dydx[2703] = 1;
}