#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_kowald1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_0000001;
    y[1] = species_0000002;
    y[2] = species_0000006;
    y[3] = species_0000007;
    y[4] = species_0000008;
    y[5] = species_0000009;
    y[6] = species_0000011;
    y[7] = species_0000016;
    y[8] = species_0000017;
}