#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_bindschadler1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.29999999999999999;
    x0[1] = 0.10000000000000001;
    x0[2] = 0.80000000000000004;
    x0[3] = 0.10000000000000001;
}