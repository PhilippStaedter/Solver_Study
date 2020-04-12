#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_kolodkin2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Lc;
    y[1] = Ln_;
    y[2] = NRLn;
    y[3] = NRn;
    y[4] = RE;
    y[5] = REL;
}