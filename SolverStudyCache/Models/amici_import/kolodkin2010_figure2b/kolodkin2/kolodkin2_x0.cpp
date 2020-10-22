#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_kolodkin2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.0050000000000000001;
    x0[3] = 370.0;
    x0[4] = 3.7000000000000002;
}