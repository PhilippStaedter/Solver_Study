#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Bungay2006a(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1174.5;
    x0[4] = 1.75;
    x0[8] = 7.5999999999999996;
    x0[10] = 0.10000000000000001;
    x0[12] = 142.84999999999999;
    x0[18] = 116.0;
    x0[22] = 66.0;
    x0[24] = 0.018200000000000001;
    x0[37] = 12.300000000000001;
    x0[38] = 4721.0;
    x0[46] = 1.0;
    x0[50] = 170000.0;
    x0[51] = 364.0;
}