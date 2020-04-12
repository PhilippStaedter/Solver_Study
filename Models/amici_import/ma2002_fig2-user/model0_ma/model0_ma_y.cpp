#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_ma(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ACA;
    y[1] = CAR1;
    y[2] = ERK2;
    y[3] = PKA;
    y[4] = REGA;
    y[5] = excAMP;
    y[6] = incAMP;
}