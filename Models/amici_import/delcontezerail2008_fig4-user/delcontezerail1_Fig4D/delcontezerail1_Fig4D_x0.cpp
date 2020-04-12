#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_delcontezerail1_Fig4D(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.001;
    x0[1] = 0.001;
    x0[2] = 1.0;
    x0[3] = 1.0;
}