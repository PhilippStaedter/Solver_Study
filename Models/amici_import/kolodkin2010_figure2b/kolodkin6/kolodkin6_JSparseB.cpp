#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_kolodkin6(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = -1.0*dwdx0;
    JSparseB[1] = -1.0*dwdx2 + 1.0*dwdx3;
    JSparseB[2] = 1.0*dwdx4;
    JSparseB[3] = 1.0*dwdx9;
    JSparseB[4] = 1.0*dwdx3;
    JSparseB[5] = 1.0*dwdx4 - 3.4444444444444402*dwdx5;
    JSparseB[6] = 1.0*dwdx9;
    JSparseB[7] = -3.4444444444444402*dwdx12;
    JSparseB[8] = -1.0*dwdx1;
    JSparseB[9] = -1.0*dwdx6 + 1.0*dwdx7;
    JSparseB[10] = 1.0*dwdx10;
    JSparseB[11] = -1.0*dwdx11;
    JSparseB[12] = -1.0*dwdx3;
    JSparseB[13] = -1.0*dwdx4;
    JSparseB[14] = -3.4444444444444402*dwdx7;
    JSparseB[15] = -3.4444444444444402*dwdx10 + 1.0*dwdx8 - 1.0*dwdx9;
    JSparseB[16] = 1.0*dwdx13;
    JSparseB[17] = 1.0*dwdx14;
    JSparseB[18] = 1.0*dwdx1;
    JSparseB[19] = 1.0*dwdx5;
    JSparseB[20] = 1.0*dwdx6;
    JSparseB[21] = 1.0*dwdx11 + 1.0*dwdx12;
    JSparseB[22] = 1.0*dwdx8;
    JSparseB[23] = 1.0*dwdx13;
    JSparseB[24] = 1.0*dwdx14;
    JSparseB[25] = -1.0*dwdx8;
    JSparseB[26] = -1.0*dwdx13;
    JSparseB[27] = -1.0*dwdx14;
}