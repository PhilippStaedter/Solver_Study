#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_gorlich1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.6000000000000001;
    x0[1] = 500.0;
    x0[2] = 0.69999999999999996;
    x0[6] = 2.0;
    x0[7] = 0.69999999999999996;
    x0[8] = 5.0;
}