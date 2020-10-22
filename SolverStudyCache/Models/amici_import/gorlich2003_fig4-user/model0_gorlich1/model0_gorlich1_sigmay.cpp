#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void sigmay_model0_gorlich1(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
    sigmay[1] = 1.0;
    sigmay[2] = 1.0;
    sigmay[3] = 1.0;
    sigmay[4] = 1.0;
    sigmay[5] = 1.0;
    sigmay[6] = 1.0;
    sigmay[7] = 1.0;
    sigmay[8] = 1.0;
    sigmay[9] = 1.0;
    sigmay[10] = 1.0;
    sigmay[11] = 1.0;
    sigmay[12] = 1.0;
}