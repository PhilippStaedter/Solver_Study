#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_bachmann2(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 2.5*dwdx1;
    JSparseB[1] = -2.5*dwdx3;
    JSparseB[2] = 2.5*dwdx2;
    JSparseB[3] = -2.5*dwdx8;
    JSparseB[4] = 3.6363636363636362*dwdx4;
    JSparseB[5] = -3.6363636363636362*dwdx43;
    JSparseB[6] = -3.6363636363636362*dwdx4;
    JSparseB[7] = 3.6363636363636362*dwdx5;
    JSparseB[8] = -3.6363636363636362*dwdx5;
    JSparseB[9] = 3.6363636363636362*dwdx6;
    JSparseB[10] = -3.6363636363636362*dwdx6;
    JSparseB[11] = 3.6363636363636362*dwdx7;
    JSparseB[12] = -3.6363636363636362*dwdx7;
    JSparseB[13] = 3.6363636363636362*dwdx8;
    JSparseB[14] = 2.5*dwdx9;
    JSparseB[15] = 2.5*dwdx10;
    JSparseB[16] = -2.5*dwdx16;
    JSparseB[17] = -2.5*dwdx21 - 2.5*dwdx22 - 2.5*dwdx23 - 2.5*dwdx24;
    JSparseB[18] = 2.5*dwdx25;
    JSparseB[19] = -2.5*dwdx49;
    JSparseB[20] = -2.5*dwdx55;
    JSparseB[21] = -2.5*dwdx59;
    JSparseB[22] = 2.5*dwdx11;
    JSparseB[23] = 2.5*dwdx45;
    JSparseB[24] = 2.5*dwdx50;
    JSparseB[25] = -2.5*dwdx9;
    JSparseB[26] = -2.5*dwdx10;
    JSparseB[27] = 2.5*dwdx12;
    JSparseB[28] = 2.5*dwdx16 + 2.5*dwdx17 + 2.5*dwdx18;
    JSparseB[29] = 2.5*dwdx21;
    JSparseB[30] = -2.5*dwdx25 + 2.5*dwdx28 + 2.5*dwdx30;
    JSparseB[31] = 2.5*dwdx14;
    JSparseB[32] = 2.5*dwdx19;
    JSparseB[33] = -2.5*dwdx20;
    JSparseB[34] = 2.5*dwdx46;
    JSparseB[35] = 2.5*dwdx51;
    JSparseB[36] = 2.5*dwdx56;
    JSparseB[37] = -2.5*dwdx14;
    JSparseB[38] = -2.5*dwdx19;
    JSparseB[39] = 2.5*dwdx20;
    JSparseB[40] = -2.5*dwdx46;
    JSparseB[41] = -2.5*dwdx51;
    JSparseB[42] = -2.5*dwdx56;
    JSparseB[43] = 2.5*dwdx29;
    JSparseB[44] = -2.5*dwdx34;
    JSparseB[45] = 2.5*dwdx33;
    JSparseB[46] = -2.5*dwdx39;
    JSparseB[47] = 3.6363636363636362*dwdx35;
    JSparseB[48] = -3.6363636363636362*dwdx44;
    JSparseB[49] = -3.6363636363636362*dwdx35;
    JSparseB[50] = 3.6363636363636362*dwdx36;
    JSparseB[51] = -3.6363636363636362*dwdx36;
    JSparseB[52] = 3.6363636363636362*dwdx37;
    JSparseB[53] = -3.6363636363636362*dwdx37;
    JSparseB[54] = 3.6363636363636362*dwdx38;
    JSparseB[55] = -3.6363636363636362*dwdx38;
    JSparseB[56] = 3.6363636363636362*dwdx39;
    JSparseB[57] = 2.5*dwdx0;
    JSparseB[58] = 2.5*dwdx15;
    JSparseB[59] = 2.5*dwdx26 + 2.5*dwdx27;
    JSparseB[60] = 2.5*dwdx40 + 2.5*dwdx41;
    JSparseB[61] = -2.5*dwdx42;
    JSparseB[62] = 2.5*dwdx47 + 2.5*dwdx48;
    JSparseB[63] = 2.5*dwdx52 + 2.5*dwdx53;
    JSparseB[64] = 2.5*dwdx57;
    JSparseB[65] = 3.6363636363636362*dwdx42;
    JSparseB[66] = -3.6363636363636362*dwdx60;
    JSparseB[67] = -2.5*dwdx13;
    JSparseB[68] = 2.5*dwdx24;
    JSparseB[69] = -2.5*dwdx31 - 2.5*dwdx32;
    JSparseB[70] = 2.5*dwdx49;
    JSparseB[71] = -2.5*dwdx54;
    JSparseB[72] = -2.5*dwdx58;
    JSparseB[73] = 2.5*dwdx13;
    JSparseB[74] = -2.5*dwdx17;
    JSparseB[75] = 2.5*dwdx22;
    JSparseB[76] = -2.5*dwdx28 + 2.5*dwdx31;
    JSparseB[77] = 2.5*dwdx54 + 2.5*dwdx55;
    JSparseB[78] = -2.5*dwdx12;
    JSparseB[79] = -2.5*dwdx18;
    JSparseB[80] = 2.5*dwdx23;
    JSparseB[81] = -2.5*dwdx30 + 2.5*dwdx32;
    JSparseB[82] = 2.5*dwdx58 + 2.5*dwdx59;
    JSparseB[83] = -2.5*dwdx0;
    JSparseB[84] = -2.5*dwdx15;
    JSparseB[85] = -2.5*dwdx26 - 2.5*dwdx27;
    JSparseB[86] = -2.5*dwdx40 - 2.5*dwdx41;
    JSparseB[87] = -2.5*dwdx47 - 2.5*dwdx48;
    JSparseB[88] = -2.5*dwdx52 - 2.5*dwdx53;
    JSparseB[89] = -2.5*dwdx57;
    JSparseB[90] = 2.5*dwdx60;
}