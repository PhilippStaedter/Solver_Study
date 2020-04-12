#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Ueda2001(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 0.29999999999999999;
    x0[2] = 0.40000000000000002;
    x0[3] = 0.20000000000000001;
    x0[4] = 0.10000000000000001;
    x0[5] = 0.59999999999999998;
    x0[6] = 0.5;
    x0[7] = 0.90000000000000002;
    x0[8] = 1.0;
    x0[9] = 0.80000000000000004;
    x0[10] = 0.69999999999999996;
    x0[11] = 1.0;
    x0[12] = 1.0;
}