#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model1_kouril4(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ADP;
    y[1] = ATP;
    y[2] = BPG;
    y[3] = P3G;
    y[4] = pep;
    y[5] = pyr;
}