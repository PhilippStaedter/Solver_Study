#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_kolodkin3(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.0050000000000000001;
    x0[5] = 42.0;
    x0[6] = 42.0;
    x0[7] = 83.25;
    x0[8] = 3.7000000000000002;
}