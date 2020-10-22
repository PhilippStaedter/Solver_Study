#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Sasagawa2005(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[100] = 1;
    dydx[200] = 1;
    dydx[300] = 1;
    dydx[400] = 1;
    dydx[500] = 1;
    dydx[600] = 1;
    dydx[700] = 1;
    dydx[800] = 1;
    dydx[900] = 1;
    dydx[1000] = 1;
    dydx[1100] = 1;
    dydx[1200] = 1;
    dydx[1300] = 1;
    dydx[1400] = 1;
    dydx[1500] = 1;
    dydx[1600] = 1;
    dydx[1700] = 1;
    dydx[1800] = 1;
    dydx[1900] = 1;
    dydx[2000] = 1;
    dydx[2100] = 1;
    dydx[2200] = 1;
    dydx[2300] = 1;
    dydx[2400] = 1;
    dydx[2500] = 1;
    dydx[2600] = 1;
    dydx[2700] = 1;
    dydx[2800] = 1;
    dydx[2900] = 1;
    dydx[3000] = 1;
    dydx[3100] = 1;
    dydx[3200] = 1;
    dydx[3300] = 1;
    dydx[3400] = 1;
    dydx[3500] = 1;
    dydx[3600] = 1;
    dydx[3700] = 1;
    dydx[3800] = 1;
    dydx[3900] = 1;
    dydx[4000] = 1;
    dydx[4100] = 1;
    dydx[4200] = 1;
    dydx[4300] = 1;
    dydx[4400] = 1;
    dydx[4500] = 1;
    dydx[4600] = 1;
    dydx[4700] = 1;
    dydx[4800] = 1;
    dydx[4900] = 1;
    dydx[5000] = 1;
    dydx[5100] = 1;
    dydx[5200] = 1;
    dydx[5300] = 1;
    dydx[5400] = 1;
    dydx[5500] = 1;
    dydx[5600] = 1;
    dydx[5700] = 1;
    dydx[5800] = 1;
    dydx[5900] = 1;
    dydx[6000] = 1;
    dydx[6100] = 1;
    dydx[6200] = 1;
    dydx[6300] = 1;
    dydx[6400] = 1;
    dydx[6500] = 1;
    dydx[6600] = 1;
    dydx[6700] = 1;
    dydx[6800] = 1;
    dydx[6900] = 1;
    dydx[7000] = 1;
    dydx[7100] = 1;
    dydx[7200] = 1;
    dydx[7300] = 1;
    dydx[7400] = 1;
    dydx[7500] = 1;
    dydx[7600] = 1;
    dydx[7700] = 1;
    dydx[7800] = 1;
    dydx[7900] = 1;
    dydx[8000] = 1;
    dydx[8100] = 1;
    dydx[8200] = 1;
    dydx[8300] = 1;
    dydx[8400] = 1;
    dydx[8500] = 1;
    dydx[8600] = 1;
    dydx[8700] = 1;
    dydx[8800] = 1;
    dydx[8900] = 1;
    dydx[9000] = 1;
    dydx[9100] = 1;
    dydx[9200] = 1;
    dydx[9300] = 1;
    dydx[9400] = 1;
    dydx[9500] = 1;
    dydx[9600] = 1;
    dydx[9700] = 1;
    dydx[9800] = 1;
}