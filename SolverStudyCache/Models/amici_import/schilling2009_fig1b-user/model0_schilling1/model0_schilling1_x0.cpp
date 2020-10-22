#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_schilling1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[8] = 7.0;
    x0[9] = 21.0;
    x0[10] = 50.0;
    x0[11] = 1.0;
    x0[12] = 2.0;
    x0[13] = 24.0;
    x0[14] = 11.0;
    x0[15] = 3.7719;
    x0[16] = 10.799099999999999;
    x0[17] = 2.5101;
}