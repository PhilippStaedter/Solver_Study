#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_BIOMD0000000037(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 10.0;
    x0[5] = 30.0;
    x0[7] = 6.0;
    x0[8] = 0.90000000000000002;
    x0[11] = 200.0;
}