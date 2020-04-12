#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_bachmann(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[7] = 1.2499699999999999e-7;
    x0[8] = 3.9762200000000001;
    x0[11] = 26.725100000000001;
    x0[20] = 79.753500000000003;
}