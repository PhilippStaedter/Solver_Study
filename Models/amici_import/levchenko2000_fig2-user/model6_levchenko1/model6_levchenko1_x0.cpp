#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model6_levchenko1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.40000000000000002;
    x0[2] = 0.29999999999999999;
    x0[8] = 0.20000000000000001;
    x0[9] = 0.20000000000000001;
    x0[16] = 0.29999999999999999;
    x0[17] = 0.20000000000000001;
    x0[18] = 0.29999999999999999;
}