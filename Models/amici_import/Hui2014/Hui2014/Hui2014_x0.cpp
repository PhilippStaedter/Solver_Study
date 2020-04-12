#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Hui2014(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[3] = 500.0;
    x0[5] = 500.0;
    x0[8] = 90.0;
    x0[9] = 10.0;
    x0[11] = 30.0;
    x0[12] = 25.0;
    x0[14] = 75.0;
    x0[17] = 100.0;
    x0[21] = 100.0;
    x0[23] = 40.0;
    x0[24] = 360.0;
    x0[27] = 1500.0;
    x0[30] = 100.0;
    x0[35] = 2.0;
    x0[37] = 100.0;
    x0[38] = 600.0;
    x0[41] = 600.0;
    x0[44] = 600.0;
    x0[46] = 2.0;
    x0[47] = 100.0;
    x0[49] = 10.0;
    x0[57] = 1000.0;
    x0[61] = 200.0;
    x0[63] = 1.0;
}