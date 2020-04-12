#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_xie1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2;
    JSparseB[1] = -1.0*dwdx8;
    JSparseB[2] = -1.0*dwdx9;
    JSparseB[3] = -1.0*dwdx11;
    JSparseB[4] = 1.0*dwdx16;
    JSparseB[5] = -1.0*dwdx0;
    JSparseB[6] = 1.0*dwdx7 + 1.0*dwdx8;
    JSparseB[7] = -1.0*dwdx16;
    JSparseB[8] = -1.0*dwdx2;
    JSparseB[9] = 1.0*dwdx10 + 1.0*dwdx9;
    JSparseB[10] = 1.0*dwdx11;
    JSparseB[11] = -1.0*dwdx23;
    JSparseB[12] = 1.0*dwdx12;
    JSparseB[13] = -1.0*dwdx27;
    JSparseB[14] = 1.0*dwdx14 + 1.0*dwdx15;
    JSparseB[15] = -1.0*dwdx18;
    JSparseB[16] = 1.0*dwdx19;
    JSparseB[17] = -1.0*dwdx30;
    JSparseB[18] = 1.0*dwdx0;
    JSparseB[19] = -1.0*dwdx8;
    JSparseB[20] = -1.0*dwdx14;
    JSparseB[21] = 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18;
    JSparseB[22] = -1.0*dwdx19;
    JSparseB[23] = 1.0*dwdx14;
    JSparseB[24] = -1.0*dwdx18;
    JSparseB[25] = 1.0*dwdx19 + 1.0*dwdx20;
    JSparseB[26] = -1.0*dwdx53;
    JSparseB[27] = 1.0*dwdx21;
    JSparseB[28] = -1.0*dwdx56;
    JSparseB[29] = 1.0*dwdx24;
    JSparseB[30] = -1.0*dwdx25;
    JSparseB[31] = -1.0*dwdx44;
    JSparseB[32] = -1.0*dwdx48;
    JSparseB[33] = 1.0*dwdx26;
    JSparseB[34] = -1.0*dwdx28;
    JSparseB[35] = -1.0*dwdx32;
    JSparseB[36] = 1.0*dwdx29;
    JSparseB[37] = -1.0*dwdx31;
    JSparseB[38] = -1.0*dwdx35;
    JSparseB[39] = -1.0*dwdx5;
    JSparseB[40] = 1.0*dwdx33 - 1.0*dwdx34;
    JSparseB[41] = -1.0*dwdx3;
    JSparseB[42] = -1.0*dwdx36 + 1.0*dwdx37;
    JSparseB[43] = -1.0*dwdx6;
    JSparseB[44] = -1.0*dwdx39 + 1.0*dwdx40;
    JSparseB[45] = -1.0*dwdx4;
    JSparseB[46] = 1.0*dwdx42 - 1.0*dwdx43;
    JSparseB[47] = -1.0*dwdx13;
    JSparseB[48] = -1.0*dwdx46 + 1.0*dwdx47;
    JSparseB[49] = -1.0*dwdx51;
    JSparseB[50] = -1.0*dwdx22;
    JSparseB[51] = -1.0*dwdx45;
    JSparseB[52] = -1.0*dwdx49 + 1.0*dwdx50;
    JSparseB[53] = -1.0*dwdx38;
    JSparseB[54] = 1.0*dwdx52;
    JSparseB[55] = -1.0*dwdx54;
    JSparseB[56] = -1.0*dwdx41;
    JSparseB[57] = 1.0*dwdx55;
    JSparseB[58] = -1.0*dwdx57;
}