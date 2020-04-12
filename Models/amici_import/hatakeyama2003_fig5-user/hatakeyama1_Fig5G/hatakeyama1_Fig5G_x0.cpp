#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_hatakeyama1_Fig5G(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 10.0;
    x0[4] = 7.0;
    x0[5] = 1000.0;
    x0[8] = 10.0;
    x0[9] = 330.0;
    x0[10] = 120.0;
    x0[13] = 2.3999999999999999;
    x0[14] = 10.0;
    x0[17] = 11.4;
    x0[18] = 800.0;
    x0[19] = 80.0;
    x0[28] = 100.0;
    x0[30] = 120.0;
    x0[34] = 1000.0;
}