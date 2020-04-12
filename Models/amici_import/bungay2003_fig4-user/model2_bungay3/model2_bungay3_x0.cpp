#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model2_bungay3(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[5] = 3400.0;
    x0[6] = 1400.0;
    x0[12] = 90.0;
    x0[19] = 849079.0;
    x0[20] = 60.0;
    x0[22] = 300.0;
    x0[26] = 2.5;
    x0[34] = 0.0050000000000000001;
    x0[35] = 1.0;
    x0[38] = 0.69999999999999996;
    x0[46] = 10.0;
    x0[48] = 0.10000000000000001;
    x0[52] = 20.0;
    x0[60] = 30.0;
    x0[63] = 170.0;
}