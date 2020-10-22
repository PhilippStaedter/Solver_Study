#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_kholodenko2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 280.0;
    x0[1] = 10.0;
    x0[2] = 10.0;
    x0[3] = 280.0;
    x0[4] = 90.0;
    x0[5] = 10.0;
    x0[6] = 10.0;
    x0[7] = 10.0;
}