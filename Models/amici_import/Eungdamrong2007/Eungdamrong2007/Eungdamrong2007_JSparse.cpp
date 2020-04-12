#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_Eungdamrong2007(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 20.833333333333332*dwdx0 - 20.833333333333332*dwdx1 - 20.833333333333332*dwdx2 - 20.833333333333332*dwdx3;
    JSparse[1] = 1.0245901639344261*dwdx1;
    JSparse[2] = 20.833333333333332*dwdx2;
    JSparse[3] = -1.0245901639344261*dwdx0;
    JSparse[4] = -1.0245901639344261*dwdx3;
    JSparse[5] = 20.833333333333332*dwdx3;
    JSparse[6] = -1.6666666666666667*dwdx4;
    JSparse[7] = 1.6666666666666667*dwdx4;
    JSparse[8] = -1.0245901639344261*dwdx5;
    JSparse[9] = -1.0245901639344261*dwdx5;
    JSparse[10] = 1.0245901639344261*dwdx5;
    JSparse[11] = -1.0245901639344261*dwdx6;
    JSparse[12] = 83.333333333333329*dwdx6;
    JSparse[13] = -1.6666666666666667*dwdx7 + 1.6666666666666667*dwdx8;
    JSparse[14] = -1.6666666666666667*dwdx8;
    JSparse[15] = 1.6666666666666667*dwdx7;
    JSparse[16] = 1.6666666666666667*dwdx9;
    JSparse[17] = -1.6666666666666667*dwdx9;
    JSparse[18] = -1.6666666666666667*dwdx10 + 1.6666666666666667*dwdx11;
    JSparse[19] = 1.6666666666666667*dwdx10 - 1.6666666666666667*dwdx11;
    JSparse[20] = 1.6666666666666667*dwdx12;
    JSparse[21] = -1.0245901639344261*dwdx12;
    JSparse[22] = -1.6666666666666667*dwdx13;
    JSparse[23] = 1.6666666666666667*dwdx13;
    JSparse[24] = 21.929824561403507*dwdx14;
    JSparse[25] = -1.0245901639344261*dwdx14 - 1.0245901639344261*dwdx15;
    JSparse[26] = -21.929824561403507*dwdx14;
    JSparse[27] = 83.333333333333329*dwdx15;
    JSparse[28] = 1.6666666666666667*dwdx17;
    JSparse[29] = -1.6666666666666667*dwdx16 - 1.6666666666666667*dwdx17;
    JSparse[30] = -1.0245901639344261*dwdx16;
    JSparse[31] = 1.6666666666666667*dwdx16;
    JSparse[32] = -1.6666666666666667*dwdx18;
    JSparse[33] = 1.6666666666666667*dwdx18;
    JSparse[34] = -1.6666666666666667*dwdx19;
    JSparse[35] = 1.6666666666666667*dwdx19;
    JSparse[36] = -1.6666666666666667*dwdx20;
    JSparse[37] = 1.0245901639344261*dwdx20;
    JSparse[38] = 1.6666666666666667*dwdx20;
    JSparse[39] = -20.833333333333332*dwdx22;
    JSparse[40] = -1.0245901639344261*dwdx21 + 1.0245901639344261*dwdx22 - 1.0245901639344261*dwdx23;
    JSparse[41] = 1.0245901639344261*dwdx23;
    JSparse[42] = 1.6666666666666667*dwdx21;
    JSparse[43] = -1.6666666666666667*dwdx24 + 1.6666666666666667*dwdx25;
    JSparse[44] = 1.6666666666666667*dwdx24 - 1.6666666666666667*dwdx25;
    JSparse[45] = -1.6666666666666667*dwdx26;
    JSparse[46] = 1.6666666666666667*dwdx26 - 1.6666666666666667*dwdx27;
    JSparse[47] = 1.0245901639344261*dwdx27;
    JSparse[48] = 1.6666666666666667*dwdx27;
    JSparse[49] = -1.6666666666666667*dwdx28;
    JSparse[50] = 1.6666666666666667*dwdx28;
    JSparse[51] = -1.6666666666666667*dwdx30;
    JSparse[52] = 1.6666666666666667*dwdx31;
    JSparse[53] = -1.6666666666666667*dwdx31;
    JSparse[54] = -1.6666666666666667*dwdx29 + 1.6666666666666667*dwdx30;
    JSparse[55] = -1.0245901639344261*dwdx32;
    JSparse[56] = 21.929824561403507*dwdx37;
    JSparse[57] = -1.0245901639344261*dwdx32 - 1.0245901639344261*dwdx33 - 1.0245901639344261*dwdx34 - 1.0245901639344261*dwdx35 - 1.0245901639344261*dwdx36 - 1.0245901639344261*dwdx37 - 1.0245901639344261*dwdx38 - 1.0245901639344261*dwdx39 - 1.0245901639344261*dwdx40;
    JSparse[58] = -21.929824561403507*dwdx37;
    JSparse[59] = -21.929824561403507*dwdx36;
    JSparse[60] = 21.929824561403507*dwdx36;
    JSparse[61] = 1.0245901639344261*dwdx32;
    JSparse[62] = 1.0245901639344261*dwdx35;
    JSparse[63] = 83.333333333333329*dwdx38 + 83.333333333333329*dwdx39 + 83.333333333333329*dwdx40;
    JSparse[64] = -1.0245901639344261*dwdx35;
    JSparse[65] = -1.0245901639344261*dwdx34;
    JSparse[66] = 1.0245901639344261*dwdx34;
    JSparse[67] = -1.0245901639344261*dwdx33;
    JSparse[68] = 1.0245901639344261*dwdx33;
    JSparse[69] = 21.929824561403507*dwdx41;
    JSparse[70] = -1.0245901639344261*dwdx41 - 1.0245901639344261*dwdx42;
    JSparse[71] = -21.929824561403507*dwdx41;
    JSparse[72] = 83.333333333333329*dwdx42;
    JSparse[73] = -1.0245901639344261*dwdx43 - 1.0245901639344261*dwdx44;
    JSparse[74] = -21.929824561403507*dwdx43;
    JSparse[75] = 21.929824561403507*dwdx43;
    JSparse[76] = 83.333333333333329*dwdx44;
    JSparse[77] = -1.0245901639344261*dwdx45 - 1.0245901639344261*dwdx46;
    JSparse[78] = -21.929824561403507*dwdx45;
    JSparse[79] = 21.929824561403507*dwdx45;
    JSparse[80] = 83.333333333333329*dwdx46;
    JSparse[81] = -1.0245901639344261*dwdx48;
    JSparse[82] = -1.0245901639344261*dwdx47;
    JSparse[83] = 83.333333333333329*dwdx48;
    JSparse[84] = 20.833333333333332*dwdx49 + 20.833333333333332*dwdx52;
    JSparse[85] = -20.833333333333332*dwdx49 - 20.833333333333332*dwdx50 + 20.833333333333332*dwdx51 - 20.833333333333332*dwdx52;
    JSparse[86] = -1.0245901639344261*dwdx51;
    JSparse[87] = 1.0245901639344261*dwdx50;
    JSparse[88] = 20.833333333333332*dwdx53;
    JSparse[89] = -20.833333333333332*dwdx53;
    JSparse[90] = 20.833333333333332*dwdx54;
    JSparse[91] = -1.0245901639344261*dwdx54;
    JSparse[92] = -20.833333333333332*dwdx55 - 20.833333333333332*dwdx56;
    JSparse[93] = 20.833333333333332*dwdx55;
    JSparse[94] = -1.0245901639344261*dwdx55;
    JSparse[95] = 20.833333333333332*dwdx57;
    JSparse[96] = -20.833333333333332*dwdx57;
    JSparse[97] = -20.833333333333332*dwdx58;
    JSparse[98] = 20.833333333333332*dwdx58;
    JSparse[99] = -1.0245901639344261*dwdx58;
    JSparse[100] = -1.0245901639344261*dwdx59;
    JSparse[101] = 1.6666666666666667*dwdx60;
    JSparse[102] = -1.0245901639344261*dwdx59;
    JSparse[103] = 1.0245901639344261*dwdx59 - 1.0245901639344261*dwdx60;
    JSparse[104] = 20.833333333333332*dwdx61;
    JSparse[105] = -1.0245901639344261*dwdx61 - 1.0245901639344261*dwdx62 + 1.0245901639344261*dwdx63;
    JSparse[106] = 1.0245901639344261*dwdx62;
    JSparse[107] = -1.6666666666666667*dwdx63;
    JSparse[108] = 20.833333333333332*dwdx64;
    JSparse[109] = -1.0245901639344261*dwdx65;
    JSparse[110] = -1.0245901639344261*dwdx64 + 1.0245901639344261*dwdx65 + 1.0245901639344261*dwdx66;
    JSparse[111] = -1.6666666666666667*dwdx66;
    JSparse[112] = -1.0245901639344261*dwdx69;
    JSparse[113] = -20.833333333333332*dwdx67;
    JSparse[114] = 1.0245901639344261*dwdx67 - 1.0245901639344261*dwdx68 + 1.0245901639344261*dwdx69;
    JSparse[115] = 1.6666666666666667*dwdx68;
    JSparse[116] = -20.833333333333332*dwdx71;
    JSparse[117] = -1.0245901639344261*dwdx70;
    JSparse[118] = 1.0245901639344261*dwdx70 - 1.0245901639344261*dwdx71;
    JSparse[119] = 20.833333333333332*dwdx71;
    JSparse[120] = -1.0245901639344261*dwdx70;
    JSparse[121] = -20.833333333333332*dwdx73;
    JSparse[122] = 1.0245901639344261*dwdx72;
    JSparse[123] = 20.833333333333332*dwdx72;
    JSparse[124] = -1.0245901639344261*dwdx73;
    JSparse[125] = 20.833333333333332*dwdx73;
    JSparse[126] = 1.0245901639344261*dwdx74;
    JSparse[127] = 20.833333333333332*dwdx74;
    JSparse[128] = -1.0245901639344261*dwdx75 - 1.0245901639344261*dwdx76 - 1.0245901639344261*dwdx77;
    JSparse[129] = 83.333333333333329*dwdx75 + 83.333333333333329*dwdx76 + 83.333333333333329*dwdx77;
    JSparse[130] = -1.0245901639344261*dwdx78 - 1.0245901639344261*dwdx79;
    JSparse[131] = 83.333333333333329*dwdx78 + 83.333333333333329*dwdx79;
    JSparse[132] = -1.0245901639344261*dwdx80;
    JSparse[133] = -1.0245901639344261*dwdx80;
    JSparse[134] = 1.0245901639344261*dwdx80;
    JSparse[135] = -1.0245901639344261*dwdx81;
    JSparse[136] = -1.0245901639344261*dwdx81;
    JSparse[137] = 1.0245901639344261*dwdx81;
    JSparse[138] = -1.0245901639344261*dwdx82;
    JSparse[139] = 1.0245901639344261*dwdx82;
    JSparse[140] = -1.0245901639344261*dwdx82;
    JSparse[141] = -1.0245901639344261*dwdx83;
    JSparse[142] = -1.0245901639344261*dwdx83;
    JSparse[143] = 1.0245901639344261*dwdx83;
    JSparse[144] = -1.0245901639344261*dwdx84;
    JSparse[145] = -1.0245901639344261*dwdx84;
    JSparse[146] = 1.0245901639344261*dwdx84;
    JSparse[147] = -1.6666666666666667*dwdx86;
    JSparse[148] = -1.0245901639344261*dwdx85;
    JSparse[149] = -1.0245901639344261*dwdx85;
    JSparse[150] = 1.0245901639344261*dwdx85 - 1.0245901639344261*dwdx86;
    JSparse[151] = 1.6666666666666667*dwdx86;
    JSparse[152] = -1.6666666666666667*dwdx87;
    JSparse[153] = -1.0245901639344261*dwdx87;
    JSparse[154] = 1.6666666666666667*dwdx87;
    JSparse[155] = 1.6666666666666667*dwdx88;
    JSparse[156] = -1.6666666666666667*dwdx88;
    JSparse[157] = -1.0245901639344261*dwdx89;
    JSparse[158] = 1.0245901639344261*dwdx92;
    JSparse[159] = 1.6666666666666667*dwdx89 - 1.6666666666666667*dwdx90 - 1.6666666666666667*dwdx91 - 1.6666666666666667*dwdx92;
    JSparse[160] = 1.6666666666666667*dwdx90 + 1.6666666666666667*dwdx91;
    JSparse[161] = 1.0245901639344261*dwdx95;
    JSparse[162] = -1.0245901639344261*dwdx94;
    JSparse[163] = -1.6666666666666667*dwdx93 + 1.6666666666666667*dwdx96;
    JSparse[164] = 1.6666666666666667*dwdx93 + 1.6666666666666667*dwdx94 - 1.6666666666666667*dwdx95 - 1.6666666666666667*dwdx96;
    JSparse[165] = -1.0245901639344261*dwdx97;
    JSparse[166] = -20.833333333333332*dwdx98;
    JSparse[167] = 20.833333333333332*dwdx98;
    JSparse[168] = -1.0245901639344261*dwdx97 - 1.0245901639344261*dwdx98;
    JSparse[169] = 1.0245901639344261*dwdx97;
    JSparse[170] = -1.0245901639344261*dwdx99;
    JSparse[171] = 20.833333333333332*dwdx100;
    JSparse[172] = -1.0245901639344261*dwdx99;
    JSparse[173] = -1.0245901639344261*dwdx100 + 1.0245901639344261*dwdx99;
}