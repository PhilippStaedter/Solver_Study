#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_marhl(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.76000000000000001;
    x0[1] = 0.28999999999999998;
    x0[2] = 85.450000000000003;
    x0[3] = 0.34999999999999998;
    x0[4] = 34.549999999999997;
}