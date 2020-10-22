#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_tripathi1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 500.0;
    x0[1] = 1.0;
    x0[2] = 10000.0;
    x0[3] = 2000.0;
    x0[4] = 2500.0;
}