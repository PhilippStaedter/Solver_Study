#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_miao2008(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 32554830.0;
    x0[1] = 134173.0;
    x0[2] = 9818.0;
    x0[3] = 26180.0;
}