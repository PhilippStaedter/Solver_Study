#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_fisher1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 9.0999999999999993e-6;
    x0[1] = 5.0500000000000001e-5;
    x0[2] = 1.0;
    x0[3] = 1.0;
    x0[4] = 0.0097108000000000003;
    x0[5] = 0.049197999999999999;
    x0[6] = 6.1e-6;
    x0[7] = 0.0009477;
    x0[8] = 0.00011010000000000001;
    x0[9] = 0.00052189999999999995;
    x0[10] = 2.2000000000000001e-6;
    x0[11] = 2.5000000000000002e-6;
    x0[12] = 0.0094397000000000005;
    x0[13] = 0.00022719999999999999;
}