#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model3_mellor1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_1;
    y[1] = species_10;
    y[2] = species_11;
    y[3] = species_12;
    y[4] = species_13;
    y[5] = species_14;
    y[6] = species_15;
    y[7] = species_7;
    y[8] = species_8;
    y[9] = species_9;
}