#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_lee2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = E;
    y[1] = E_M;
    y[2] = E_P2;
    y[3] = E_P_1;
    y[4] = E_P_2;
    y[5] = M;
    y[6] = P;
    y[7] = P2;
    y[8] = T;
}