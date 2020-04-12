#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_saeidi1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.0000000000000001e-5;
    x0[1] = 1.0000000000000001e-5;
    x0[3] = 9.9999999999999995e-7;
    x0[6] = 5.0000000000000004e-6;
}