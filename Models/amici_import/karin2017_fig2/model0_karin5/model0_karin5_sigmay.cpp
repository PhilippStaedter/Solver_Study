#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void sigmay_model0_karin5(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
    sigmay[1] = 1.0;
    sigmay[2] = 1.0;
    sigmay[3] = 1.0;
    sigmay[4] = 1.0;
}