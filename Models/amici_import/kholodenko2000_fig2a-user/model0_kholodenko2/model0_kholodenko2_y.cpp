#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_kholodenko2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = MAPK;
    y[1] = MAPK_P;
    y[2] = MAPK_PP;
    y[3] = MKK;
    y[4] = MKKK;
    y[5] = MKKK_P;
    y[6] = MKK_P;
    y[7] = MKK_PP;
}