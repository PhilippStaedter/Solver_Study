#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Bungay2003(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1400.0;
    x0[4] = 20.0;
    x0[8] = 10.0;
    x0[10] = 0.10000000000000001;
    x0[12] = 0.69999999999999996;
    x0[16] = 90.0;
    x0[20] = 170.0;
    x0[26] = 300.0;
    x0[32] = 60.0;
    x0[34] = 0.0050000000000000001;
    x0[52] = 30.0;
    x0[57] = 2.5;
    x0[58] = 3400.0;
    x0[68] = 1.0;
    x0[73] = 849079.0;
}