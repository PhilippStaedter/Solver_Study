#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_ODea2007(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[4] = 0.059999999999999998;
    x0[8] = 0.001;
    x0[13] = 0.029999999999999999;
    x0[20] = 0.01;
}