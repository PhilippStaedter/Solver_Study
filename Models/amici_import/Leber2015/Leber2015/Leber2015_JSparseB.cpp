#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_Leber2015(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1 - 1.0*dwdx3 + 1.0*dwdx7;
    JSparseB[1] = 1.0*dwdx10;
    JSparseB[2] = 1.0*dwdx13;
    JSparseB[3] = 1.0*dwdx19 - 1.0*dwdx20;
    JSparseB[4] = 1.0*dwdx24;
    JSparseB[5] = 1.0*dwdx29;
    JSparseB[6] = 1.0*dwdx35;
    JSparseB[7] = 1.0*dwdx39;
    JSparseB[8] = 1.0*dwdx9;
    JSparseB[9] = 1.0*dwdx12;
    JSparseB[10] = 1.0*dwdx26;
    JSparseB[11] = 1.0*dwdx34;
    JSparseB[12] = -1.0*dwdx9;
    JSparseB[13] = -1.0*dwdx12 + 1.0*dwdx15;
    JSparseB[14] = -1.0*dwdx26;
    JSparseB[15] = -1.0*dwdx34;
    JSparseB[16] = -1.0*dwdx7;
    JSparseB[17] = -1.0*dwdx10;
    JSparseB[18] = -1.0*dwdx13;
    JSparseB[19] = 1.0*dwdx16;
    JSparseB[20] = -1.0*dwdx29;
    JSparseB[21] = -1.0*dwdx35;
    JSparseB[22] = -1.0*dwdx16;
    JSparseB[23] = 1.0*dwdx17 + 1.0*dwdx18;
    JSparseB[24] = 1.0*dwdx21;
    JSparseB[25] = 1.0*dwdx36;
    JSparseB[26] = 1.0*dwdx45;
    JSparseB[27] = -1.0*dwdx2;
    JSparseB[28] = 1.0*dwdx8;
    JSparseB[29] = 1.0*dwdx22;
    JSparseB[30] = -1.0*dwdx31;
    JSparseB[31] = -1.0*dwdx44;
    JSparseB[32] = -1.0*dwdx48;
    JSparseB[33] = -1.0*dwdx54;
    JSparseB[34] = 0.25*dwdx5;
    JSparseB[35] = 0.25*dwdx23;
    JSparseB[36] = 0.25*dwdx27 + 0.25*dwdx28;
    JSparseB[37] = -0.25*dwdx32;
    JSparseB[38] = 0.25*dwdx38;
    JSparseB[39] = 0.25*dwdx47;
    JSparseB[40] = -0.25*dwdx23 - 0.25*dwdx25;
    JSparseB[41] = -0.25*dwdx27;
    JSparseB[42] = 0.25*dwdx32;
    JSparseB[43] = -0.25*dwdx33 - 0.25*dwdx37;
    JSparseB[44] = -0.25*dwdx38 - 0.25*dwdx40;
    JSparseB[45] = -0.25*dwdx47 - 0.25*dwdx50;
    JSparseB[46] = -0.25*dwdx5;
    JSparseB[47] = 0.25*dwdx25;
    JSparseB[48] = -0.25*dwdx28;
    JSparseB[49] = 0.25*dwdx33 + 0.25*dwdx37;
    JSparseB[50] = 0.25*dwdx40;
    JSparseB[51] = 0.25*dwdx50;
    JSparseB[52] = -0.25*dwdx6;
    JSparseB[53] = 0.25*dwdx41;
    JSparseB[54] = -0.25*dwdx43;
    JSparseB[55] = -0.25*dwdx51;
    JSparseB[56] = -0.25*dwdx56;
    JSparseB[57] = -14.285714285714285*dwdx0;
    JSparseB[58] = 14.285714285714285*dwdx42;
    JSparseB[59] = 14.285714285714285*dwdx4;
    JSparseB[60] = 14.285714285714285*dwdx46 + 14.285714285714285*dwdx49;
    JSparseB[61] = 14.285714285714285*dwdx55;
    JSparseB[62] = -14.285714285714285*dwdx61;
    JSparseB[63] = 14.285714285714285*dwdx52;
    JSparseB[64] = -14.285714285714285*dwdx62;
    JSparseB[65] = -14.285714285714285*dwdx4;
    JSparseB[66] = -14.285714285714285*dwdx49;
    JSparseB[67] = 14.285714285714285*dwdx53 - 14.285714285714285*dwdx55;
    JSparseB[68] = -14.285714285714285*dwdx60;
    JSparseB[69] = 1.0*dwdx11;
    JSparseB[70] = 1.0*dwdx14;
    JSparseB[71] = 1.0*dwdx30;
    JSparseB[72] = -1.0*dwdx42;
    JSparseB[73] = 1.0*dwdx57 + 1.0*dwdx58 + 1.0*dwdx59;
    JSparseB[74] = -1.0*dwdx18;
    JSparseB[75] = 1.0*dwdx60;
    JSparseB[76] = -1.0*dwdx58;
    JSparseB[77] = 1.0*dwdx61;
    JSparseB[78] = -1.0*dwdx11;
    JSparseB[79] = -1.0*dwdx14;
    JSparseB[80] = -1.0*dwdx30;
    JSparseB[81] = -1.0*dwdx59;
    JSparseB[82] = 1.0*dwdx62;
}