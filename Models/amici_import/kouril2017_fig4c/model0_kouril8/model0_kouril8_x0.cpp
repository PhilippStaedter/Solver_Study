#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_kouril8(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 10.0;
    x0[4] = 3.0;
    x0[6] = 10.0;
    x0[8] = 0.20000000000000001;
    x0[9] = 5.0;
    x0[11] = 0.40000000000000002;
}