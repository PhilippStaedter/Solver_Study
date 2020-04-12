#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model3_sarma7(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_0;
    y[1] = species_1;
    y[2] = species_10;
    y[3] = species_11;
    y[4] = species_12;
    y[5] = species_2;
    y[6] = species_3;
    y[7] = species_4;
    y[8] = species_5;
    y[9] = species_6;
    y[10] = species_7;
    y[11] = species_8;
    y[12] = species_9;
}