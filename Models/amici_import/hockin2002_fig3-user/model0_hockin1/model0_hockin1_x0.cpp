#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_hockin1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 3.4000000000000001e-6;
    x0[1] = 1.3999999999999999e-6;
    x0[4] = 8.9999999999999999e-8;
    x0[9] = 2.5000000000000001e-11;
    x0[10] = 2.5000000000000001e-9;
    x0[18] = 2.0e-8;
    x0[19] = 1.0e-8;
    x0[20] = 6.9999999999999996e-10;
    x0[24] = 1.0e-10;
    x0[26] = 1.6e-7;
}