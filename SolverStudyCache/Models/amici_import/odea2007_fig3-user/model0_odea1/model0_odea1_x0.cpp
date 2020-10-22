#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_odea1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.001;
    x0[3] = 0.059999999999999998;
    x0[10] = 0.029999999999999999;
    x0[17] = 0.01;
}