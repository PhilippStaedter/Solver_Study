#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_becker1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Epo;
    y[1] = EpoR;
    y[2] = Epo_EpoR;
    y[3] = Epo_EpoRi;
    y[4] = dEpoe;
    y[5] = dEpoi;
}