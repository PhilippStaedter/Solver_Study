#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_laub1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 1.0;
    x0[2] = 1.0;
    x0[3] = 2.5;
    x0[4] = 1.3999999999999999;
    x0[5] = 1.5;
    x0[6] = 1.6000000000000001;
}