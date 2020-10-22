#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Yang2007(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = x1;
    y[1] = x10;
    y[2] = x11;
    y[3] = x12;
    y[4] = x13;
    y[5] = x14;
    y[6] = x15;
    y[7] = x16;
    y[8] = x17;
    y[9] = x18;
    y[10] = x19;
    y[11] = x2;
    y[12] = x20;
    y[13] = x21;
    y[14] = x22;
    y[15] = x23;
    y[16] = x24;
    y[17] = x25;
    y[18] = x3;
    y[19] = x4;
    y[20] = x5;
    y[21] = x6;
    y[22] = x7;
    y[23] = x8;
    y[24] = x9;
}