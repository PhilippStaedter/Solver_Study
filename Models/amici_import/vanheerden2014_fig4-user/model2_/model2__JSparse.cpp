#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model2_(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -2.0*dwdx0 + 4.0*dwdx1 - 1.0*dwdx2;
    JSparse[1] = 1.0*dwdx0 - 1.0*dwdx1;
    JSparse[2] = -2.0*dwdx1 + 1.0*dwdx2;
    JSparse[3] = 4.0*dwdx3;
    JSparse[4] = -1.0*dwdx3 - 1.0*dwdx4;
    JSparse[5] = -2.0*dwdx3 + 2.0*dwdx4;
    JSparse[6] = 4.0*dwdx5;
    JSparse[7] = -1.0*dwdx5;
    JSparse[8] = -2.0*dwdx5 + 1.0*dwdx6;
}