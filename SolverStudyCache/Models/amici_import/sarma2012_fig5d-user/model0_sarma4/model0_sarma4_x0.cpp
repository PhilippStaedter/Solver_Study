#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_sarma4(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 300.0;
    x0[1] = 200.0;
    x0[2] = 20.0;
    x0[4] = 1199.9999422132501;
    x0[7] = 1199.9999422132501;
    x0[10] = 100.0;
}