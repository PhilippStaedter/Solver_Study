#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_chassagnole2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.16700000000000001;
    x0[1] = 0.098000000000000004;
    x0[2] = 0.59999999999999998;
    x0[3] = 0.27200000000000002;
    x0[4] = 0.65300000000000002;
    x0[5] = 3.48;
    x0[6] = 0.218;
    x0[7] = 2.0;
    x0[8] = 2.6699999999999999;
    x0[9] = 0.80800000000000005;
    x0[10] = 0.39900000000000002;
    x0[11] = 2.1299999999999999;
    x0[12] = 0.0080000000000000002;
    x0[13] = 2.6699999999999999;
    x0[14] = 0.39800000000000002;
    x0[15] = 0.111;
    x0[16] = 0.27600000000000002;
    x0[17] = 0.13800000000000001;
}