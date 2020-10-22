#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_bungay1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[6] = 3400.0;
    x0[7] = 1400.0;
    x0[13] = 90.0;
    x0[20] = 170000.0;
    x0[21] = 60.0;
    x0[23] = 300.0;
    x0[27] = 2.5;
    x0[35] = 0.0050000000000000001;
    x0[36] = 1.0;
    x0[39] = 0.69999999999999996;
    x0[47] = 10.0;
    x0[49] = 0.10000000000000001;
    x0[53] = 20.0;
    x0[61] = 30.0;
    x0[64] = 170.0;
    x0[74] = 2600.0;
}