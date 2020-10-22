#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Lai2014(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[144] = 3.3000000000000003e-5;
    x0[193] = 0.000146;
}