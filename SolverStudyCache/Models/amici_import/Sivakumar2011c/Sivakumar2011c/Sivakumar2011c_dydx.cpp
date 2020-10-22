#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Sivakumar2011c(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[51] = 1;
    dydx[102] = 1;
    dydx[153] = 1;
    dydx[204] = 1;
    dydx[255] = 1;
    dydx[306] = 1;
    dydx[357] = 1;
    dydx[408] = 1;
    dydx[459] = 1;
    dydx[510] = 1;
    dydx[561] = 1;
    dydx[612] = 1;
    dydx[663] = 1;
    dydx[714] = 1;
    dydx[765] = 1;
    dydx[816] = 1;
    dydx[867] = 1;
    dydx[918] = 1;
    dydx[969] = 1;
    dydx[1020] = 1;
    dydx[1071] = 1;
    dydx[1122] = 1;
    dydx[1173] = 1;
    dydx[1224] = 1;
    dydx[1275] = 1;
    dydx[1326] = 1;
    dydx[1377] = 1;
    dydx[1428] = 1;
    dydx[1479] = 1;
    dydx[1530] = 1;
    dydx[1581] = 1;
    dydx[1632] = 1;
    dydx[1683] = 1;
    dydx[1734] = 1;
    dydx[1785] = 1;
    dydx[1836] = 1;
    dydx[1887] = 1;
    dydx[1938] = 1;
    dydx[1989] = 1;
    dydx[2040] = 1;
    dydx[2091] = 1;
    dydx[2142] = 1;
    dydx[2193] = 1;
    dydx[2244] = 1;
    dydx[2295] = 1;
    dydx[2346] = 1;
    dydx[2397] = 1;
    dydx[2448] = 1;
    dydx[2499] = 1;
}