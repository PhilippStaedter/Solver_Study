#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_piedrafita1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = S;
    y[1] = ST;
    y[2] = STU;
    y[3] = STUS;
    y[4] = STUST;
    y[5] = STUSU;
    y[6] = SU;
    y[7] = SUST;
    y[8] = SUSTU;
    y[9] = T;
    y[10] = U;
}