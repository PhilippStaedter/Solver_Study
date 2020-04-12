#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_leloup2_Fig2A(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = species_0;
    y[1] = species_1;
    y[2] = species_10;
    y[3] = species_11;
    y[4] = species_12;
    y[5] = species_13;
    y[6] = species_14;
    y[7] = species_15;
    y[8] = species_2;
    y[9] = species_3;
    y[10] = species_4;
    y[11] = species_5;
    y[12] = species_6;
    y[13] = species_7;
    y[14] = species_8;
    y[15] = species_9;
}