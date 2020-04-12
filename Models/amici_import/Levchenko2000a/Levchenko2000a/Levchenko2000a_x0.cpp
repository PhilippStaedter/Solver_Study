#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Levchenko2000a(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.29999999999999999;
    x0[1] = 0.20000000000000001;
    x0[2] = 0.10000000000000001;
    x0[3] = 0.29999999999999999;
    x0[4] = 0.40000000000000002;
    x0[7] = 0.20000000000000001;
    x0[10] = 0.29999999999999999;
    x0[22] = 0.10000000000000001;
}