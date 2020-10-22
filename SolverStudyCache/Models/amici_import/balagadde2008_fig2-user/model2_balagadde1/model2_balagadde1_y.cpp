#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model2_balagadde1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = A1;
    y[1] = A2;
    y[2] = C1;
    y[3] = C2;
    y[4] = IPTG;
    y[5] = sink;
    y[6] = source;
}