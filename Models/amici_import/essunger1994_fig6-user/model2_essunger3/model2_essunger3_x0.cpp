#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model2_essunger3(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 10.0;
    x0[2] = 540.0;
    x0[4] = 1000.0;
    x0[5] = 450.0;
    x0[6] = 5.0000000000000002e-5;
}