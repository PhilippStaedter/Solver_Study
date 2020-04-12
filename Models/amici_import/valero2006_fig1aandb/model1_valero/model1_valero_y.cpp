#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model1_valero(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ADP;
    y[1] = AMP;
    y[2] = ATP;
    y[3] = Lac;
    y[4] = NADH;
    y[5] = Pyr;
}