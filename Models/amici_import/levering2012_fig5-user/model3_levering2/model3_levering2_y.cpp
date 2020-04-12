#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model3_levering2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_1;
    y[1] = species_10;
    y[2] = species_11;
    y[3] = species_12;
    y[4] = species_13;
    y[5] = species_14;
    y[6] = species_15;
    y[7] = species_16;
    y[8] = species_17;
    y[9] = species_18;
    y[10] = species_19;
    y[11] = species_2;
    y[12] = species_20;
    y[13] = species_21;
    y[14] = species_22;
    y[15] = species_23;
    y[16] = species_24;
    y[17] = species_25;
    y[18] = species_3;
    y[19] = species_4;
    y[20] = species_5;
    y[21] = species_6;
    y[22] = species_7;
    y[23] = species_8;
    y[24] = species_9;
}