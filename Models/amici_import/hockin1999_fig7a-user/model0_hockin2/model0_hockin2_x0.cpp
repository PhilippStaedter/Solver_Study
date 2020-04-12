#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_hockin2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.0e-8;
    x0[10] = 1.9999999999999999e-7;
}