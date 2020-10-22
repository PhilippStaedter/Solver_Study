#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_neves1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 300.0;
    x0[3] = 3010000.0;
    x0[4] = 3010000.0;
    x0[6] = 94.0;
    x0[8] = 120.40000000000001;
    x0[10] = 0.60199999999999998;
    x0[13] = 2167.1999999999998;
    x0[15] = 216.72;
    x0[17] = 108.36;
    x0[19] = 240.80000000000001;
    x0[20] = 301.0;
    x0[22] = 60.200000000000003;
    x0[23] = 120.40000000000001;
    x0[25] = 60.200000000000003;
    x0[26] = 120.40000000000001;
    x0[27] = 120.40000000000001;
    x0[36] = 6020.0;
}