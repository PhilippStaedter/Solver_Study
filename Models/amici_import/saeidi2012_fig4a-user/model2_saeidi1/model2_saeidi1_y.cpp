#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model2_saeidi1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = s1;
    y[1] = s16;
    y[2] = s17;
    y[3] = s19;
    y[4] = s2;
    y[5] = s3;
    y[6] = s4;
    y[7] = s42;
    y[8] = s45;
    y[9] = s5;
}