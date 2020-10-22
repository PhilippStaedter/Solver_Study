#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_hornberg1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Inh;
    y[1] = R;
    y[2] = Rin;
    y[3] = x1;
    y[4] = x1p;
    y[5] = x2;
    y[6] = x2p;
    y[7] = x3;
    y[8] = x3p;
}