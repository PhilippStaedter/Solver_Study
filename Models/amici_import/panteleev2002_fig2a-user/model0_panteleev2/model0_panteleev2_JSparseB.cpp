#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_panteleev2(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JSparseB[1] = 1.0*dwdx9;
    JSparseB[2] = 1.0*dwdx10;
    JSparseB[3] = 1.0*dwdx15;
    JSparseB[4] = 1.0*dwdx16;
    JSparseB[5] = 1.0*dwdx2 - 1.0*dwdx3 + 1.0*dwdx4 - 1.0*dwdx5;
    JSparseB[6] = 1.0*dwdx6;
    JSparseB[7] = -1.0*dwdx8;
    JSparseB[8] = -1.0*dwdx11;
    JSparseB[9] = 1.0*dwdx13;
    JSparseB[10] = -1.0*dwdx14;
    JSparseB[11] = 1.0*dwdx17 - 1.0*dwdx18;
    JSparseB[12] = 1.0*dwdx19;
    JSparseB[13] = -1.0*dwdx2;
    JSparseB[14] = -1.0*dwdx6 + 1.0*dwdx7;
    JSparseB[15] = -1.0*dwdx13;
    JSparseB[16] = 1.0*dwdx1;
    JSparseB[17] = 1.0*dwdx3;
    JSparseB[18] = -1.0*dwdx7;
    JSparseB[19] = 1.0*dwdx8 + 1.0*dwdx9;
    JSparseB[20] = 1.0*dwdx10;
    JSparseB[21] = 1.0*dwdx14;
    JSparseB[22] = -1.0*dwdx1;
    JSparseB[23] = 1.0*dwdx5;
    JSparseB[24] = -1.0*dwdx9;
    JSparseB[25] = -1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12;
    JSparseB[26] = 1.0*dwdx18;
    JSparseB[27] = 1.0*dwdx20;
    JSparseB[28] = 1.0*dwdx2;
    JSparseB[29] = 1.0*dwdx6;
    JSparseB[30] = 1.0*dwdx13;
    JSparseB[31] = 1.0*dwdx0;
    JSparseB[32] = -1.0*dwdx3;
    JSparseB[33] = -1.0*dwdx8;
    JSparseB[34] = -1.0*dwdx14 + 1.0*dwdx15;
    JSparseB[35] = 1.0*dwdx16;
    JSparseB[36] = -1.0*dwdx0;
    JSparseB[37] = 1.0*dwdx4 - 1.0*dwdx5;
    JSparseB[38] = -1.0*dwdx11;
    JSparseB[39] = -1.0*dwdx15;
    JSparseB[40] = -1.0*dwdx16 + 1.0*dwdx17 - 1.0*dwdx18;
    JSparseB[41] = 1.0*dwdx19;
    JSparseB[42] = -1.0*dwdx4;
    JSparseB[43] = -1.0*dwdx12;
    JSparseB[44] = -1.0*dwdx17;
    JSparseB[45] = -1.0*dwdx19 - 1.0*dwdx20;
}