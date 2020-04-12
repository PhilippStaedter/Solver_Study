#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_bungay2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[4] = 4721.0;
    x0[5] = 1174.5;
    x0[11] = 170000.0;
    x0[12] = 66.0;
    x0[14] = 116.0;
    x0[18] = 12.300000000000001;
    x0[24] = 0.018200000000000001;
    x0[25] = 1.0;
    x0[27] = 7.5999999999999996;
    x0[29] = 0.10000000000000001;
    x0[33] = 1.75;
    x0[40] = 142.84999999999999;
    x0[50] = 364.0;
}