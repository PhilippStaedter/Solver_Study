#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_ihekwaba1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    JSparse[1] = 1.0*dwdx5;
    JSparse[2] = 1.0*dwdx2;
    JSparse[3] = 1.0*dwdx0;
    JSparse[4] = 1.0*dwdx3;
    JSparse[5] = 1.0*dwdx1;
    JSparse[6] = 1.0*dwdx4;
    JSparse[7] = -1.0*dwdx5;
    JSparse[8] = -1.0*dwdx2;
    JSparse[9] = -1.0*dwdx0;
    JSparse[10] = -1.0*dwdx3;
    JSparse[11] = -1.0*dwdx1;
    JSparse[12] = -1.0*dwdx4;
    JSparse[13] = 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[14] = -1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JSparse[15] = 1.0*dwdx7;
    JSparse[16] = 1.0*dwdx9;
    JSparse[17] = -1.0*dwdx7;
    JSparse[18] = 1.0*dwdx10 + 1.0*dwdx11;
    JSparse[19] = 1.0*dwdx12;
    JSparse[20] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    JSparse[21] = 1.0*dwdx11;
    JSparse[22] = 1.0*dwdx10 + 1.0*dwdx12;
    JSparse[23] = 1.0*dwdx13 + 1.0*dwdx14;
    JSparse[24] = -1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15;
    JSparse[25] = 1.0*dwdx15;
    JSparse[26] = 1.0*dwdx14;
    JSparse[27] = -1.0*dwdx15;
    JSparse[28] = 1.0*dwdx16 + 1.0*dwdx17;
    JSparse[29] = 1.0*dwdx18;
    JSparse[30] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    JSparse[31] = 1.0*dwdx16;
    JSparse[32] = 1.0*dwdx17 + 1.0*dwdx18;
    JSparse[33] = 1.0*dwdx20 + 1.0*dwdx21;
    JSparse[34] = -1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21;
    JSparse[35] = 1.0*dwdx19;
    JSparse[36] = 1.0*dwdx21;
    JSparse[37] = -1.0*dwdx19;
    JSparse[38] = 1.0*dwdx22 + 1.0*dwdx24;
    JSparse[39] = 1.0*dwdx23;
    JSparse[40] = -1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JSparse[41] = 1.0*dwdx24;
    JSparse[42] = 1.0*dwdx22 + 1.0*dwdx23;
    JSparse[43] = -1.0*dwdx27;
    JSparse[44] = 1.0*dwdx27;
    JSparse[45] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JSparse[46] = 1.0*dwdx25;
    JSparse[47] = 1.0*dwdx26;
    JSparse[48] = -1.0*dwdx25;
    JSparse[49] = -1.0*dwdx29;
    JSparse[50] = 1.0*dwdx29;
    JSparse[51] = 1.0*dwdx30;
    JSparse[52] = -1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    JSparse[53] = 1.0*dwdx30 + 1.0*dwdx31;
    JSparse[54] = 1.0*dwdx33;
    JSparse[55] = -1.0*dwdx32 - 1.0*dwdx33;
    JSparse[56] = 1.0*dwdx32;
    JSparse[57] = -1.0*dwdx32;
    JSparse[58] = 1.0*dwdx35;
    JSparse[59] = 1.0*dwdx34;
    JSparse[60] = -1.0*dwdx34 - 1.0*dwdx35;
    JSparse[61] = 1.0*dwdx34;
    JSparse[62] = 1.0*dwdx37;
    JSparse[63] = -1.0*dwdx36;
    JSparse[64] = -1.0*dwdx38;
    JSparse[65] = 1.0*dwdx38;
    JSparse[66] = -1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41;
    JSparse[67] = 1.0*dwdx39;
    JSparse[68] = 1.0*dwdx40;
    JSparse[69] = -1.0*dwdx39;
    JSparse[70] = -1.0*dwdx42;
    JSparse[71] = 1.0*dwdx42;
    JSparse[72] = 1.0*dwdx43;
    JSparse[73] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    JSparse[74] = 1.0*dwdx43 + 1.0*dwdx44;
    JSparse[75] = 1.0*dwdx46;
    JSparse[76] = -1.0*dwdx45 - 1.0*dwdx46;
    JSparse[77] = 1.0*dwdx45;
    JSparse[78] = -1.0*dwdx45;
    JSparse[79] = 1.0*dwdx48;
    JSparse[80] = 1.0*dwdx47;
    JSparse[81] = -1.0*dwdx47 - 1.0*dwdx48;
    JSparse[82] = 1.0*dwdx47;
    JSparse[83] = 1.0*dwdx50;
    JSparse[84] = -1.0*dwdx49;
    JSparse[85] = -1.0*dwdx51;
    JSparse[86] = 1.0*dwdx51;
    JSparse[87] = -1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JSparse[88] = 1.0*dwdx52;
    JSparse[89] = 1.0*dwdx53;
    JSparse[90] = -1.0*dwdx52;
    JSparse[91] = -1.0*dwdx55;
    JSparse[92] = 1.0*dwdx55;
    JSparse[93] = 1.0*dwdx56;
    JSparse[94] = -1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57;
    JSparse[95] = 1.0*dwdx56 + 1.0*dwdx57;
    JSparse[96] = 1.0*dwdx59;
    JSparse[97] = -1.0*dwdx58 - 1.0*dwdx59;
    JSparse[98] = 1.0*dwdx58;
    JSparse[99] = -1.0*dwdx58;
    JSparse[100] = 1.0*dwdx61;
    JSparse[101] = 1.0*dwdx60;
    JSparse[102] = -1.0*dwdx60 - 1.0*dwdx61;
    JSparse[103] = 1.0*dwdx60;
    JSparse[104] = 1.0*dwdx63;
    JSparse[105] = -1.0*dwdx62;
    JSparse[106] = -1.0*dwdx64;
    JSparse[107] = 1.0*dwdx64;
    JSparse[108] = -1.0*dwdx70;
    JSparse[109] = 1.0*dwdx70;
    JSparse[110] = -1.0*dwdx65;
    JSparse[111] = 1.0*dwdx65;
    JSparse[112] = -1.0*dwdx66;
    JSparse[113] = 1.0*dwdx66;
    JSparse[114] = -1.0*dwdx67;
    JSparse[115] = 1.0*dwdx67;
    JSparse[116] = -1.0*dwdx68;
    JSparse[117] = 1.0*dwdx68;
    JSparse[118] = -1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66 - 1.0*dwdx67 - 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70;
    JSparse[119] = 1.0*dwdx69;
    JSparse[120] = -1.0*dwdx72;
    JSparse[121] = 1.0*dwdx72;
    JSparse[122] = 1.0*dwdx71;
    JSparse[123] = -1.0*dwdx73;
    JSparse[124] = 1.0*dwdx73;
    JSparse[125] = -1.0*dwdx74;
    JSparse[126] = 1.0*dwdx74;
    JSparse[127] = 1.0*dwdx75;
    JSparse[128] = -1.0*dwdx72 - 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75;
    JSparse[129] = 1.0*dwdx76;
    JSparse[130] = 1.0*dwdx77;
    JSparse[131] = 1.0*dwdx78;
}