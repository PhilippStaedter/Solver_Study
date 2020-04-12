#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_ferreira1_Fig2F(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[2] = 0.25;
    x0[4] = 0.0033999999999999998;
}