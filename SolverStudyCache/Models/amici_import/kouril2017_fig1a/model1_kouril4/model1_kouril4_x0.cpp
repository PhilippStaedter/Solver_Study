#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_kouril4(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 10.0;
    x0[3] = 4.0;
    x0[4] = 5.0;
}