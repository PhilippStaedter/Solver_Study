#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_karin5(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 400.0;
    x0[1] = 4.9666670000000002;
    x0[2] = 11.42;
    x0[4] = 0.27000000000000002;
}