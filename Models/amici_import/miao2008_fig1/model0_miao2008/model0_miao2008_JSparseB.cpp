#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_miao2008(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = -1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JSparseB[1] = 1.0*dwdx4;
    JSparseB[2] = 1.0*dwdx8;
    JSparseB[3] = 1.0*dwdx10;
    JSparseB[4] = -1.0*dwdx1 - 0.25*dwdx3;
    JSparseB[5] = -1.0*dwdx4 - 1.0*dwdx5 + 1.0*dwdx6;
    JSparseB[6] = -0.25*dwdx8;
    JSparseB[7] = 1.0*dwdx11;
    JSparseB[8] = -0.5*dwdx3;
    JSparseB[9] = -1.0*dwdx6 - 1.0*dwdx7;
    JSparseB[10] = -0.5*dwdx8 - 1.0*dwdx9;
    JSparseB[11] = -1.0*dwdx11 - 1.0*dwdx13;
    JSparseB[12] = -1.0*dwdx2 - 0.25*dwdx3;
    JSparseB[13] = 1.0*dwdx7;
    JSparseB[14] = -0.25*dwdx8;
    JSparseB[15] = -1.0*dwdx10 - 1.0*dwdx12 + 1.0*dwdx13;
}