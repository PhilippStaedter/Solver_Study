#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Leber2015(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 484.0;
    x0[1] = 1.0;
    x0[2] = 50000000000.0;
    x0[5] = 15000000000.0;
    x0[7] = 1052500.0;
    x0[9] = 500000.0;
    x0[11] = 3250.0;
    x0[13] = 1714285.7142857099;
    x0[14] = 714285.71428571397;
    x0[20] = 12000000.0;
}