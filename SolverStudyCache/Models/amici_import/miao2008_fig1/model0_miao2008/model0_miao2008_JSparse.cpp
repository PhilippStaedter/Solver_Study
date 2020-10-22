#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_miao2008(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparse[1] = 1.0*dwdx1 + 0.25*dwdx3;
    JSparse[2] = 0.5*dwdx3;
    JSparse[3] = 1.0*dwdx2 + 0.25*dwdx3;
    JSparse[4] = -1.0*dwdx4;
    JSparse[5] = 1.0*dwdx4 + 1.0*dwdx5 - 1.0*dwdx6;
    JSparse[6] = 1.0*dwdx6 + 1.0*dwdx7;
    JSparse[7] = -1.0*dwdx7;
    JSparse[8] = -1.0*dwdx8;
    JSparse[9] = 0.25*dwdx8;
    JSparse[10] = 0.5*dwdx8 + 1.0*dwdx9;
    JSparse[11] = 0.25*dwdx8;
    JSparse[12] = -1.0*dwdx10;
    JSparse[13] = -1.0*dwdx11;
    JSparse[14] = 1.0*dwdx11 + 1.0*dwdx13;
    JSparse[15] = 1.0*dwdx10 + 1.0*dwdx12 - 1.0*dwdx13;
}