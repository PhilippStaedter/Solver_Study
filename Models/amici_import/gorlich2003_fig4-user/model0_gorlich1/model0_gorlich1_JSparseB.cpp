#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_gorlich1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 83333333333.333328*dwdx2 - 83333333333.333328*dwdx3;
    JSparseB[1] = 83333333333.333328*dwdx7;
    JSparseB[2] = -83333333333.333328*dwdx9;
    JSparseB[3] = 83333333333.333328*dwdx15;
    JSparseB[4] = -83333333333.333328*dwdx22;
    JSparseB[5] = -83333333333.333328*dwdx0;
    JSparseB[6] = 83333333333.333328*dwdx1;
    JSparseB[7] = -83333333333.333328*dwdx4 + 83333333333.333328*dwdx5;
    JSparseB[8] = -83333333333.333328*dwdx6;
    JSparseB[9] = 83333333333.333328*dwdx8;
    JSparseB[10] = 83333333333.333328*dwdx0;
    JSparseB[11] = -83333333333.333328*dwdx2;
    JSparseB[12] = 83333333333.333328*dwdx4;
    JSparseB[13] = 83333333333.333328*dwdx6 - 83333333333.333328*dwdx7;
    JSparseB[14] = -83333333333.333328*dwdx15;
    JSparseB[15] = -83333333333.333328*dwdx1;
    JSparseB[16] = 83333333333.333328*dwdx3;
    JSparseB[17] = -83333333333.333328*dwdx5;
    JSparseB[18] = -83333333333.333328*dwdx8 + 83333333333.333328*dwdx9;
    JSparseB[19] = 83333333333.333328*dwdx22;
    JSparseB[20] = 55555555555.555557*dwdx10;
    JSparseB[21] = -55555555555.555557*dwdx11;
    JSparseB[22] = -55555555555.555557*dwdx16 + 55555555555.555557*dwdx17;
    JSparseB[23] = 55555555555.555557*dwdx20;
    JSparseB[24] = -55555555555.555557*dwdx11 - 55555555555.555557*dwdx12;
    JSparseB[25] = -55555555555.555557*dwdx13;
    JSparseB[26] = -55555555555.555557*dwdx14;
    JSparseB[27] = -55555555555.555557*dwdx16;
    JSparseB[28] = -55555555555.555557*dwdx19;
    JSparseB[29] = 83333333333.333328*dwdx2;
    JSparseB[30] = 83333333333.333328*dwdx7;
    JSparseB[31] = 83333333333.333328*dwdx13;
    JSparseB[32] = 83333333333.333328*dwdx14 + 83333333333.333328*dwdx15;
    JSparseB[33] = -55555555555.555557*dwdx10;
    JSparseB[34] = 55555555555.555557*dwdx11;
    JSparseB[35] = 55555555555.555557*dwdx16 - 55555555555.555557*dwdx17;
    JSparseB[36] = -55555555555.555557*dwdx20;
    JSparseB[37] = 55555555555.555557*dwdx10;
    JSparseB[38] = 55555555555.555557*dwdx12;
    JSparseB[39] = 55555555555.555557*dwdx17;
    JSparseB[40] = -55555555555.555557*dwdx18 + 55555555555.555557*dwdx19 + 55555555555.555557*dwdx20;
    JSparseB[41] = -55555555555.555557*dwdx21;
    JSparseB[42] = -83333333333.333328*dwdx3;
    JSparseB[43] = -83333333333.333328*dwdx9;
    JSparseB[44] = 83333333333.333328*dwdx18;
    JSparseB[45] = 83333333333.333328*dwdx21 - 83333333333.333328*dwdx22;
}