#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_essunger1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Ta;
    y[1] = Tastarstar;
    y[2] = Tm;
    y[3] = Tmstarstar;
    y[4] = Ttot;
    y[5] = Tv;
    y[6] = V;
}