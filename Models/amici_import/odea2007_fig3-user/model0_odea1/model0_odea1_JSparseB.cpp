#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_odea1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JSparseB[1] = -1.0*dwdx7 + 1.0*dwdx9;
    JSparseB[2] = -1.0*dwdx10 + 1.0*dwdx12;
    JSparseB[3] = 1.0*dwdx14;
    JSparseB[4] = 1.0*dwdx20;
    JSparseB[5] = -1.0*dwdx28 + 1.0*dwdx30;
    JSparseB[6] = -1.0*dwdx31 + 1.0*dwdx33;
    JSparseB[7] = 1.0*dwdx35;
    JSparseB[8] = 1.0*dwdx41;
    JSparseB[9] = -1.0*dwdx49 + 1.0*dwdx51;
    JSparseB[10] = -1.0*dwdx52 + 1.0*dwdx54;
    JSparseB[11] = 1.0*dwdx56;
    JSparseB[12] = 1.0*dwdx62;
    JSparseB[13] = -1.0*dwdx4;
    JSparseB[14] = 1.0*dwdx7 + 1.0*dwdx8 - 1.0*dwdx9;
    JSparseB[15] = 1.0*dwdx11;
    JSparseB[16] = -1.0*dwdx20;
    JSparseB[17] = 1.0*dwdx70;
    JSparseB[18] = -1.0*dwdx1;
    JSparseB[19] = -1.0*dwdx8;
    JSparseB[20] = 1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    JSparseB[21] = -1.0*dwdx14;
    JSparseB[22] = -1.0*dwdx70;
    JSparseB[23] = 1.0*dwdx1;
    JSparseB[24] = 1.0*dwdx12;
    JSparseB[25] = 1.0*dwdx13 + 1.0*dwdx14 - 1.0*dwdx15;
    JSparseB[26] = -1.0*dwdx18;
    JSparseB[27] = -1.0*dwdx21;
    JSparseB[28] = -1.0*dwdx73;
    JSparseB[29] = 1.0*dwdx16 - 1.0*dwdx17 + 1.0*dwdx18;
    JSparseB[30] = -1.0*dwdx26;
    JSparseB[31] = -1.0*dwdx77;
    JSparseB[32] = 1.0*dwdx4;
    JSparseB[33] = 1.0*dwdx9;
    JSparseB[34] = 1.0*dwdx15;
    JSparseB[35] = 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22;
    JSparseB[36] = -1.0*dwdx24;
    JSparseB[37] = 1.0*dwdx27;
    JSparseB[38] = 1.0*dwdx73;
    JSparseB[39] = 1.0*dwdx23;
    JSparseB[40] = -1.0*dwdx80;
    JSparseB[41] = 1.0*dwdx17;
    JSparseB[42] = -1.0*dwdx22;
    JSparseB[43] = 1.0*dwdx25 + 1.0*dwdx26 - 1.0*dwdx27;
    JSparseB[44] = 1.0*dwdx77;
    JSparseB[45] = -1.0*dwdx5;
    JSparseB[46] = 1.0*dwdx28 + 1.0*dwdx29 - 1.0*dwdx30;
    JSparseB[47] = 1.0*dwdx32;
    JSparseB[48] = -1.0*dwdx41;
    JSparseB[49] = 1.0*dwdx71;
    JSparseB[50] = -1.0*dwdx2;
    JSparseB[51] = -1.0*dwdx29;
    JSparseB[52] = 1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33;
    JSparseB[53] = -1.0*dwdx35;
    JSparseB[54] = -1.0*dwdx71;
    JSparseB[55] = 1.0*dwdx2;
    JSparseB[56] = 1.0*dwdx33;
    JSparseB[57] = 1.0*dwdx34 + 1.0*dwdx35 - 1.0*dwdx36;
    JSparseB[58] = -1.0*dwdx39;
    JSparseB[59] = -1.0*dwdx42;
    JSparseB[60] = -1.0*dwdx74;
    JSparseB[61] = 1.0*dwdx37 - 1.0*dwdx38 + 1.0*dwdx39;
    JSparseB[62] = -1.0*dwdx47;
    JSparseB[63] = -1.0*dwdx78;
    JSparseB[64] = 1.0*dwdx5;
    JSparseB[65] = 1.0*dwdx30;
    JSparseB[66] = 1.0*dwdx36;
    JSparseB[67] = 1.0*dwdx40 + 1.0*dwdx41 + 1.0*dwdx42 + 1.0*dwdx43;
    JSparseB[68] = -1.0*dwdx45;
    JSparseB[69] = 1.0*dwdx48;
    JSparseB[70] = 1.0*dwdx74;
    JSparseB[71] = 1.0*dwdx44;
    JSparseB[72] = 1.0*dwdx38;
    JSparseB[73] = -1.0*dwdx43;
    JSparseB[74] = 1.0*dwdx46 + 1.0*dwdx47 - 1.0*dwdx48;
    JSparseB[75] = 1.0*dwdx78;
    JSparseB[76] = -1.0*dwdx6;
    JSparseB[77] = 1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51;
    JSparseB[78] = 1.0*dwdx53;
    JSparseB[79] = -1.0*dwdx62;
    JSparseB[80] = 1.0*dwdx72;
    JSparseB[81] = -1.0*dwdx3;
    JSparseB[82] = -1.0*dwdx50;
    JSparseB[83] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JSparseB[84] = -1.0*dwdx56;
    JSparseB[85] = -1.0*dwdx72;
    JSparseB[86] = 1.0*dwdx3;
    JSparseB[87] = 1.0*dwdx54;
    JSparseB[88] = 1.0*dwdx55 + 1.0*dwdx56 - 1.0*dwdx57;
    JSparseB[89] = -1.0*dwdx60;
    JSparseB[90] = -1.0*dwdx63;
    JSparseB[91] = -1.0*dwdx75;
    JSparseB[92] = 1.0*dwdx58 - 1.0*dwdx59 + 1.0*dwdx60;
    JSparseB[93] = -1.0*dwdx68;
    JSparseB[94] = -1.0*dwdx79;
    JSparseB[95] = 1.0*dwdx6;
    JSparseB[96] = 1.0*dwdx51;
    JSparseB[97] = 1.0*dwdx57;
    JSparseB[98] = 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63 + 1.0*dwdx64;
    JSparseB[99] = -1.0*dwdx66;
    JSparseB[100] = 1.0*dwdx69;
    JSparseB[101] = 1.0*dwdx75;
    JSparseB[102] = 1.0*dwdx65;
    JSparseB[103] = 1.0*dwdx59;
    JSparseB[104] = -1.0*dwdx64;
    JSparseB[105] = 1.0*dwdx67 + 1.0*dwdx68 - 1.0*dwdx69;
    JSparseB[106] = 1.0*dwdx79;
    JSparseB[107] = 1.0*dwdx8;
    JSparseB[108] = -1.0*dwdx10 + 1.0*dwdx11;
    JSparseB[109] = -1.0*dwdx13 + 1.0*dwdx15;
    JSparseB[110] = 1.0*dwdx21;
    JSparseB[111] = 1.0*dwdx29;
    JSparseB[112] = -1.0*dwdx31 + 1.0*dwdx32;
    JSparseB[113] = -1.0*dwdx34 + 1.0*dwdx36;
    JSparseB[114] = 1.0*dwdx42;
    JSparseB[115] = 1.0*dwdx50;
    JSparseB[116] = -1.0*dwdx52 + 1.0*dwdx53;
    JSparseB[117] = -1.0*dwdx55 + 1.0*dwdx57;
    JSparseB[118] = 1.0*dwdx63;
    JSparseB[119] = 1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72 + 1.0*dwdx73 + 1.0*dwdx74 + 1.0*dwdx75 + 1.0*dwdx76;
    JSparseB[120] = 1.0*dwdx81;
    JSparseB[121] = -1.0*dwdx16 + 1.0*dwdx17;
    JSparseB[122] = 1.0*dwdx26;
    JSparseB[123] = -1.0*dwdx37 + 1.0*dwdx38;
    JSparseB[124] = 1.0*dwdx47;
    JSparseB[125] = -1.0*dwdx58 + 1.0*dwdx59;
    JSparseB[126] = 1.0*dwdx68;
    JSparseB[127] = -1.0*dwdx76;
    JSparseB[128] = 1.0*dwdx77 + 1.0*dwdx78 + 1.0*dwdx79 - 1.0*dwdx81;
}