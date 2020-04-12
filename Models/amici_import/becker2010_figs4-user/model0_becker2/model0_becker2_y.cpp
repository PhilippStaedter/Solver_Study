#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_becker2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EpoR;
    y[1] = SAv;
    y[2] = SAv_EpoR;
    y[3] = SAv_EpoRi;
    y[4] = dSAve;
    y[5] = dSAvi;
}