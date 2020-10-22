#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_kowald1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 5.0000000000000004e-6;
    x0[7] = 1.0000000000000001e-5;
    x0[8] = 1.0000000000000001e-5;
}