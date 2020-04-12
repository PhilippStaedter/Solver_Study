#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Yang2007(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.001;
    x0[1] = 0.001;
    x0[2] = 0.001;
    x0[3] = 0.001;
    x0[4] = 0.001;
    x0[5] = 0.001;
    x0[6] = 1.5;
    x0[7] = 1.5;
    x0[8] = 0.5;
    x0[9] = 0.20000000000000001;
    x0[10] = 0.5;
    x0[11] = 0.001;
    x0[12] = 0.5;
    x0[13] = 5.0;
    x0[14] = 0.76000000000000001;
    x0[15] = 0.070000000000000007;
    x0[16] = 0.80000000000000004;
    x0[18] = 0.001;
    x0[19] = 0.001;
    x0[20] = 0.001;
    x0[21] = 0.001;
    x0[22] = 0.001;
    x0[23] = 0.001;
    x0[24] = 0.001;
}