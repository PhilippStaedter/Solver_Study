#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_orfao1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.3999999999999999;
    x0[4] = 0.049999990000000001;
    x0[6] = 0.029999999999999999;
    x0[7] = 0.029999990000000001;
    x0[9] = 0.19999990000000001;
}