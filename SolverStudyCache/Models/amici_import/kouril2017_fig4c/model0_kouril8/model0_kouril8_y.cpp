#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_kouril8(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ADP;
    y[1] = ATP;
    y[2] = BPG;
    y[3] = GA;
    y[4] = P3G;
    y[5] = gap;
    y[6] = glc;
    y[7] = nadp;
    y[8] = nadph;
    y[9] = pep;
    y[10] = pi;
    y[11] = pyr;
    y[12] = sink;
}