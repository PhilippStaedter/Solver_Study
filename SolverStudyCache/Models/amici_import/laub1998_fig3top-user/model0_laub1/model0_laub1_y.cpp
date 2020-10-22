#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_laub1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_0;
    y[1] = species_1;
    y[2] = species_2;
    y[3] = species_3;
    y[4] = species_4;
    y[5] = species_5;
    y[6] = species_6;
}