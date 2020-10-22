#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_mcauley1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 304.0;
    x0[1] = 100.0;
    x0[2] = 57516.0;
    x0[4] = 10000.0;
    x0[5] = 100.0;
    x0[6] = 100.0;
    x0[8] = 20.0;
    x0[9] = 100.0;
    x0[10] = 600.0;
    x0[11] = 3150.0;
    x0[13] = 20.0;
    x0[14] = 100.0;
    x0[15] = 100.0;
    x0[16] = 100.0;
    x0[17] = 100.0;
    x0[18] = 575.15999999999997;
    x0[20] = 9363.0;
    x0[23] = 45.0;
    x0[24] = 100.0;
    x0[26] = 100.0;
    x0[27] = 100.0;
    x0[28] = 400.0;
    x0[29] = 467.0;
    x0[31] = 60000.0;
}