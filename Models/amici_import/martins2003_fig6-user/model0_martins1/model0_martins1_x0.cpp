#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_martins1(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 9.0;
}