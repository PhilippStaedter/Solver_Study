#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Liu2011(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 2.0;
    x0[1] = 0.032779599999999999;
    x0[3] = 770.0;
    x0[6] = 310.0;
    x0[7] = 2470.0;
    x0[12] = 4650.0;
    x0[16] = 6.7999999999999998;
    x0[20] = 20.0;
    x0[26] = 260.0;
    x0[41] = 0.00050000000000000001;
}