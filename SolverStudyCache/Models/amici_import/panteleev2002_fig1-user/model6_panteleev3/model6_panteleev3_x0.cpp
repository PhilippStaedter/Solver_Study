#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model6_panteleev3(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 2.3999990000000002;
    x0[1] = 0.99999970000000005;
    x0[4] = 169.9999;
}