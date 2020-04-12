#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_aguda1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 1.0;
    x0[5] = 5.0;
    x0[6] = 0.01;
    x0[7] = 0.01;
    x0[8] = 2.0;
    x0[9] = 0.050000000000000003;
    x0[10] = 1.95;
}