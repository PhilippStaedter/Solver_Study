#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_liu2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 2470.0;
    x0[1] = 310.0;
    x0[4] = 4650.0;
    x0[7] = 770.0;
    x0[8] = 260.0;
    x0[17] = 2.0;
    x0[28] = 20.0;
    x0[29] = 6.7999999999999998;
    x0[30] = 0.032779599999999999;
    x0[37] = 0.00050000000000000001;
}