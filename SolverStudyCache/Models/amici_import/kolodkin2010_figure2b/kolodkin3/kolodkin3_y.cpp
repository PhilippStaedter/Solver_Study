#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_kolodkin3(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Lc;
    y[1] = Ln_;
    y[2] = NRLc;
    y[3] = NRLm;
    y[4] = NRLn;
    y[5] = NRc;
    y[6] = NRm;
    y[7] = NRn;
    y[8] = RE;
    y[9] = REL;
}