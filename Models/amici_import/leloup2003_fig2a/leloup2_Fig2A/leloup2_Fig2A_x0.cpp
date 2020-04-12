#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_leloup2_Fig2A(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 9.0;
    x0[1] = 2.0;
    x0[4] = 1.0;
    x0[9] = 1.8999999999999999;
    x0[11] = 1.3999999999999999;
    x0[13] = 1.6000000000000001;
}