#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_teusink2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 4.0;
    x0[1] = 1.0;
    x0[2] = 10.0;
    x0[3] = 0.10000000000000001;
}