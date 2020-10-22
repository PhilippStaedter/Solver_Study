#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model2_alsheihk1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1800.0;
    x0[1] = 1.0;
    x0[2] = 5400.0;
    x0[3] = 4500.0;
    x0[4] = 15300.0;
}