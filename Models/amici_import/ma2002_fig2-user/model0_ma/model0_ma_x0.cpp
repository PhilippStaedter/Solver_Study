#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_ma(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 3.3900000000000001;
    x0[1] = 2.4500000000000002;
    x0[2] = 1.1299999999999999;
    x0[3] = 1.6000000000000001;
    x0[4] = 0.90000000000000002;
    x0[5] = 0.47999999999999998;
    x0[6] = 1.2;
}