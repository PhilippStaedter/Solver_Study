#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_kholodenko1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 680.0;
    x0[2] = 85.0;
    x0[3] = 105.0;
    x0[6] = 100.0;
    x0[18] = 34.0;
    x0[22] = 150.0;
}