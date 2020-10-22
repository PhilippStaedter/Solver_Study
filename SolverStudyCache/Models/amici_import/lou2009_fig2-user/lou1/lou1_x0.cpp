#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_lou1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 39455.0;
    x0[1] = 36931.0;
    x0[2] = 234800.0;
    x0[3] = 2590800.0;
    x0[4] = 2425100.0;
    x0[5] = 15418000.0;
}