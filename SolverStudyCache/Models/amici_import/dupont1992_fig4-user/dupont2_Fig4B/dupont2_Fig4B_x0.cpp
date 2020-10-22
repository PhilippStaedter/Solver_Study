#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_dupont2_Fig4B(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 1.0;
    x0[2] = 1.6000000000000001;
    x0[3] = 0.14999999999999999;
}