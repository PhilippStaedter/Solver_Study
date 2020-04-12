#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_hornberg1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 0.5;
    x0[3] = 1.0;
    x0[5] = 1.0;
    x0[7] = 1.0;
}