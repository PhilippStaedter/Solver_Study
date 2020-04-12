#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_BIOMD0000000037(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Gluc;
    y[1] = Pfr;
    y[2] = Pi;
    y[3] = Pr;
    y[4] = S;
    y[5] = V;
    y[6] = Xa;
    y[7] = Xi;
    y[8] = Ya;
    y[9] = Yi;
    y[10] = preS;
    y[11] = prepreS;
}