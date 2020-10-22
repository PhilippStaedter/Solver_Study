#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_sarma1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1200.0;
    x0[4] = 200.0;
    x0[7] = 300.0;
    x0[9] = 10.0;
    x0[12] = 100.0;
    x0[24] = 1200.0;
}