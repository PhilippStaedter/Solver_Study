#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_kolodkin6(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.0050000000000000001;
    x0[2] = 199.40299999999999;
    x0[5] = 49.850000000000001;
    x0[6] = 3.7000000000000002;
}