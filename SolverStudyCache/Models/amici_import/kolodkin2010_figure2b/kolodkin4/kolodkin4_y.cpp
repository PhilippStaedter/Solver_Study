#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_kolodkin4(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Lc;
    y[1] = Ln_;
    y[2] = NR;
    y[3] = NRLc;
    y[4] = NRLn;
    y[5] = NRc;
    y[6] = RE;
    y[7] = REL;
}