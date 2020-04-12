#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model2_(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[3] = 1.0;
    x0[4] = 2.0;
    x0[6] = 10.0;
}