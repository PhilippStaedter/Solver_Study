#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_hald(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 1.5;
    x0[3] = 0.33000000000000002;
    x0[4] = 2.1000000000000001;
    x0[15] = 24.0;
    x0[19] = 5.0;
    x0[20] = 0.75;
    x0[21] = 0.23000000000000001;
}