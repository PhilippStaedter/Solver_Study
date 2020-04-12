#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_hornberg1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx1;
    JSparseB[1] = -1.0*dwdx3;
    JSparseB[2] = -1.0*dwdx1;
    JSparseB[3] = 1.0*dwdx3;
    JSparseB[4] = 1.0*dwdx2;
    JSparseB[5] = 1.0*dwdx4;
    JSparseB[6] = -1.0*dwdx5;
    JSparseB[7] = -1.0*dwdx2;
    JSparseB[8] = -1.0*dwdx4;
    JSparseB[9] = 1.0*dwdx5;
    JSparseB[10] = 1.0*dwdx6;
    JSparseB[11] = 1.0*dwdx7;
    JSparseB[12] = -1.0*dwdx8;
    JSparseB[13] = -1.0*dwdx6;
    JSparseB[14] = -1.0*dwdx7;
    JSparseB[15] = 1.0*dwdx8;
    JSparseB[16] = -1.0*dwdx0;
    JSparseB[17] = 1.0*dwdx9;
    JSparseB[18] = 1.0*dwdx10;
    JSparseB[19] = -1.0*dwdx11;
    JSparseB[20] = 1.0*dwdx0;
    JSparseB[21] = -1.0*dwdx9;
    JSparseB[22] = -1.0*dwdx10;
    JSparseB[23] = 1.0*dwdx11;
}