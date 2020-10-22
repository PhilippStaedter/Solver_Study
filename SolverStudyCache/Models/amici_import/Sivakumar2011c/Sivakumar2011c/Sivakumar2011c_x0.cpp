#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Sivakumar2011c(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 3.0;
    x0[4] = 3.0;
    x0[5] = 1.0;
    x0[6] = 3.0;
    x0[7] = 1.0;
    x0[10] = 1.0;
    x0[11] = 3.0;
    x0[13] = 3.0;
    x0[14] = 3.0;
    x0[15] = 4.0;
    x0[18] = 2.0;
    x0[19] = 2.0;
    x0[20] = 2.0;
    x0[21] = 2.0;
    x0[23] = 4.0;
    x0[24] = 5.0;
    x0[25] = 5.0;
    x0[35] = 5.0;
    x0[36] = 1.0;
    x0[44] = 2.0;
    x0[48] = 0.5;
}