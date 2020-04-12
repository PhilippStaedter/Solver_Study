#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_tiago1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = s1;
    y[1] = s10;
    y[2] = s11;
    y[3] = s12;
    y[4] = s13;
    y[5] = s14;
    y[6] = s15;
    y[7] = s16;
    y[8] = s17;
    y[9] = s2;
    y[10] = s3;
    y[11] = s4;
    y[12] = s5;
    y[13] = s6;
    y[14] = s7;
    y[15] = s8;
    y[16] = s9;
}