#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_zhao1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.0;
    x0[1] = 149400.0;
    x0[2] = 16600.0;
    x0[3] = 667200.0;
    x0[4] = 166800.0;
}