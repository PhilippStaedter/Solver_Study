#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_bier2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 6.5999999999999996;
    x0[1] = 10.300000000000001;
    x0[2] = 7.5999999999999996;
    x0[3] = 0.40999999999999998;
}