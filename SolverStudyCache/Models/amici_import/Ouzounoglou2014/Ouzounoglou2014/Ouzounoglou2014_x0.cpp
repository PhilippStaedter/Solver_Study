#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Ouzounoglou2014(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 200.0;
    x0[52] = 8.0;
    x0[53] = 8.0;
    x0[54] = 1.0;
    x0[55] = 8.0;
    x0[57] = 953.0;
    x0[58] = 1650.0;
    x0[59] = 22.0;
    x0[60] = 8.0;
    x0[61] = 8.0;
    x0[62] = 750.0;
    x0[63] = 8.0;
    x0[64] = 8.0;
    x0[65] = 8.0;
    x0[66] = 8.0;
    x0[67] = 8.0;
    x0[68] = 8.0;
    x0[69] = 8.0;
    x0[70] = 8.0;
    x0[71] = 8.0;
    x0[73] = 1500.0;
}