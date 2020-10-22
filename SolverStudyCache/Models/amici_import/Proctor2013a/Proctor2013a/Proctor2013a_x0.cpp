#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Proctor2013a(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[3] = 100.0;
    x0[4] = 5.0;
    x0[6] = 100.0;
    x0[9] = 100.0;
    x0[11] = 100.0;
    x0[13] = 100.0;
    x0[18] = 100.0;
    x0[27] = 100.0;
    x0[29] = 20.0;
    x0[30] = 20.0;
    x0[31] = 100.0;
    x0[37] = 100000.0;
    x0[41] = 100.0;
    x0[51] = 1000.0;
    x0[52] = 200.0;
    x0[53] = 200.0;
    x0[57] = 100.0;
    x0[62] = 100.0;
    x0[72] = 2.0;
    x0[73] = 1.0;
    x0[74] = 1.0;
}