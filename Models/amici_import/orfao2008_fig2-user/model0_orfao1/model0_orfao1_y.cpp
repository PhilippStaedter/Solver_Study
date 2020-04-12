#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_orfao1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = II;
    y[1] = IIa;
    y[2] = IIa_ATIII;
    y[3] = IIa_alpha2M;
    y[4] = PL;
    y[5] = PT;
    y[6] = RVV;
    y[7] = V;
    y[8] = Va;
    y[9] = X;
    y[10] = Xa;
    y[11] = Xa_ATIII;
}