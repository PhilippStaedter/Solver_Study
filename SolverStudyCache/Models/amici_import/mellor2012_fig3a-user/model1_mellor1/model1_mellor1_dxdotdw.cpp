#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model1_mellor1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -0.001;
    dxdotdw[1] = 1.2e-5;
    dxdotdw[2] = 0.00016200000000000001;
    dxdotdw[3] = 4.0000000000000003e-5;
    dxdotdw[4] = 1.4e-5;
    dxdotdw[5] = 3.9999999999999998e-6;
    dxdotdw[6] = 0.00057399999999999997;
    dxdotdw[7] = 0.000144;
    dxdotdw[8] = 5.0000000000000002e-5;
    dxdotdw[9] = -0.001;
    dxdotdw[10] = 1.5e-5;
    dxdotdw[11] = 0.000127;
    dxdotdw[12] = 2.5999999999999998e-5;
    dxdotdw[13] = 1.8e-5;
    dxdotdw[14] = 1.5999999999999999e-5;
    dxdotdw[15] = 0.00075100000000000004;
    dxdotdw[16] = 2.3e-5;
    dxdotdw[17] = 2.5000000000000001e-5;
    dxdotdw[18] = -0.001;
    dxdotdw[19] = 0.000107;
    dxdotdw[20] = 0.00021800000000000001;
    dxdotdw[21] = 0.00021800000000000001;
    dxdotdw[22] = 9.800000000000001e-5;
    dxdotdw[23] = 9.7e-5;
    dxdotdw[24] = 6.8000000000000013e-5;
    dxdotdw[25] = 5.8999999999999998e-5;
    dxdotdw[26] = 0.00013600000000000003;
    dxdotdw[27] = 0.001;
    dxdotdw[28] = -0.001;
    dxdotdw[29] = 0.001;
    dxdotdw[30] = -0.001;
}