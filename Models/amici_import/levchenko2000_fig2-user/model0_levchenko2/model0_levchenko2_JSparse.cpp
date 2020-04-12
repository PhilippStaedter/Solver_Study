#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_levchenko2(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0;
    JSparse[1] = 1.0*dwdx1;
    JSparse[2] = 1.0*dwdx2;
    JSparse[3] = -1.0*dwdx3;
    JSparse[4] = -1.0*dwdx4;
    JSparse[5] = 1.0*dwdx5;
    JSparse[6] = -1.0*dwdx6;
    JSparse[7] = 1.0*dwdx7;
    JSparse[8] = 1.0*dwdx8;
    JSparse[9] = -1.0*dwdx9;
    JSparse[10] = 1.0*dwdx10;
    JSparse[11] = 1.0*dwdx11;
    JSparse[12] = -1.0*dwdx12;
    JSparse[13] = 1.0*dwdx13;
    JSparse[14] = -1.0*dwdx14;
    JSparse[15] = 1.0*dwdx15;
    JSparse[16] = 1.0*dwdx16;
    JSparse[17] = -1.0*dwdx17;
    JSparse[18] = 1.0*dwdx18;
    JSparse[19] = 1.0*dwdx19;
    JSparse[20] = -1.0*dwdx20;
    JSparse[21] = 1.0*dwdx21;
    JSparse[22] = -1.0*dwdx22;
    JSparse[23] = 1.0*dwdx23;
    JSparse[24] = 1.0*dwdx24;
    JSparse[25] = -1.0*dwdx25;
    JSparse[26] = 1.0*dwdx26;
    JSparse[27] = 1.0*dwdx27;
    JSparse[28] = -1.0*dwdx28;
    JSparse[29] = 1.0*dwdx29;
    JSparse[30] = 1.0*dwdx30;
    JSparse[31] = 1.0*dwdx31;
    JSparse[32] = 1.0*dwdx32;
    JSparse[33] = 1.0*dwdx33;
    JSparse[34] = -1.0*dwdx34;
    JSparse[35] = 1.0*dwdx35;
    JSparse[36] = 1.0*dwdx36;
    JSparse[37] = 1.0*dwdx37;
    JSparse[38] = 1.0*dwdx38;
    JSparse[39] = 1.0*dwdx39;
    JSparse[40] = -1.0*dwdx40;
    JSparse[41] = 1.0*dwdx41;
    JSparse[42] = 1.0*dwdx42;
    JSparse[43] = 1.0*dwdx43;
    JSparse[44] = 1.0*dwdx44;
    JSparse[45] = 1.0*dwdx45;
    JSparse[46] = -1.0*dwdx46;
    JSparse[47] = 1.0*dwdx47;
    JSparse[48] = 1.0*dwdx48;
    JSparse[49] = -1.0*dwdx50;
    JSparse[50] = -1.0*dwdx51;
    JSparse[51] = -1.0*dwdx52;
    JSparse[52] = 1.0*dwdx53;
    JSparse[53] = 1.0*dwdx54;
    JSparse[54] = 1.0*dwdx55;
    JSparse[55] = -1.0*dwdx49 - 1.0*dwdx56;
    JSparse[56] = 1.0*dwdx49;
    JSparse[57] = -1.0*dwdx49;
    JSparse[58] = 1.0*dwdx57;
    JSparse[59] = -1.0*dwdx57 - 1.0*dwdx58;
    JSparse[60] = 1.0*dwdx58;
    JSparse[61] = 1.0*dwdx57 + 1.0*dwdx58;
    JSparse[62] = -1.0*dwdx59 - 1.0*dwdx60;
    JSparse[63] = -1.0*dwdx59;
    JSparse[64] = 1.0*dwdx59;
    JSparse[65] = -1.0*dwdx60;
    JSparse[66] = 1.0*dwdx60;
    JSparse[67] = -1.0*dwdx61;
    JSparse[68] = -1.0*dwdx61 - 1.0*dwdx62;
    JSparse[69] = 1.0*dwdx61;
    JSparse[70] = 1.0*dwdx62;
    JSparse[71] = -1.0*dwdx62;
    JSparse[72] = 1.0*dwdx64;
    JSparse[73] = 1.0*dwdx63 + 1.0*dwdx64;
    JSparse[74] = 1.0*dwdx63;
    JSparse[75] = -1.0*dwdx63 - 1.0*dwdx64;
    JSparse[76] = 1.0*dwdx65;
    JSparse[77] = -1.0*dwdx65 - 1.0*dwdx66;
    JSparse[78] = 1.0*dwdx66;
    JSparse[79] = 1.0*dwdx65 + 1.0*dwdx66;
    JSparse[80] = -1.0*dwdx67;
    JSparse[81] = -1.0*dwdx67;
    JSparse[82] = 1.0*dwdx67;
    JSparse[83] = 1.0*dwdx68 + 1.0*dwdx69;
    JSparse[84] = 1.0*dwdx69;
    JSparse[85] = 1.0*dwdx68;
    JSparse[86] = -1.0*dwdx68 - 1.0*dwdx69;
    JSparse[87] = -1.0*dwdx71;
    JSparse[88] = 1.0*dwdx72;
    JSparse[89] = -1.0*dwdx73;
    JSparse[90] = -1.0*dwdx74;
    JSparse[91] = 1.0*dwdx75;
    JSparse[92] = 1.0*dwdx76;
    JSparse[93] = -1.0*dwdx70 - 1.0*dwdx77;
    JSparse[94] = 1.0*dwdx70;
    JSparse[95] = -1.0*dwdx70;
    JSparse[96] = -1.0*dwdx78 - 1.0*dwdx79;
    JSparse[97] = -1.0*dwdx78;
    JSparse[98] = 1.0*dwdx78;
    JSparse[99] = -1.0*dwdx79;
    JSparse[100] = 1.0*dwdx79;
    JSparse[101] = 1.0*dwdx80;
    JSparse[102] = -1.0*dwdx80 - 1.0*dwdx81;
    JSparse[103] = 1.0*dwdx81;
    JSparse[104] = 1.0*dwdx80 + 1.0*dwdx81;
    JSparse[105] = -1.0*dwdx82;
    JSparse[106] = -1.0*dwdx82 - 1.0*dwdx83;
    JSparse[107] = 1.0*dwdx82;
    JSparse[108] = 1.0*dwdx83;
    JSparse[109] = -1.0*dwdx83;
    JSparse[110] = 1.0*dwdx85;
    JSparse[111] = 1.0*dwdx84 + 1.0*dwdx85;
    JSparse[112] = 1.0*dwdx84;
    JSparse[113] = -1.0*dwdx84 - 1.0*dwdx85;
    JSparse[114] = 1.0*dwdx86;
    JSparse[115] = -1.0*dwdx86 - 1.0*dwdx87;
    JSparse[116] = 1.0*dwdx87;
    JSparse[117] = 1.0*dwdx86 + 1.0*dwdx87;
    JSparse[118] = -1.0*dwdx89;
    JSparse[119] = 1.0*dwdx89;
    JSparse[120] = -1.0*dwdx90;
    JSparse[121] = 1.0*dwdx90;
    JSparse[122] = -1.0*dwdx88;
    JSparse[123] = -1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90;
    JSparse[124] = 1.0*dwdx88;
    JSparse[125] = 1.0*dwdx91 + 1.0*dwdx92;
    JSparse[126] = 1.0*dwdx92;
    JSparse[127] = 1.0*dwdx91;
    JSparse[128] = -1.0*dwdx91 - 1.0*dwdx92;
    JSparse[129] = -1.0*dwdx93;
    JSparse[130] = -1.0*dwdx93;
    JSparse[131] = 1.0*dwdx93;
    JSparse[132] = -1.0*dwdx94;
    JSparse[133] = -1.0*dwdx94;
    JSparse[134] = 1.0*dwdx94;
    JSparse[135] = -1.0*dwdx95;
    JSparse[136] = -1.0*dwdx95;
    JSparse[137] = 1.0*dwdx95;
    JSparse[138] = 1.0*dwdx96;
    JSparse[139] = 1.0*dwdx96 + 1.0*dwdx97;
    JSparse[140] = -1.0*dwdx96 - 1.0*dwdx97;
    JSparse[141] = 1.0*dwdx97;
    JSparse[142] = -1.0*dwdx101;
    JSparse[143] = 1.0*dwdx102;
    JSparse[144] = -1.0*dwdx103;
    JSparse[145] = 1.0*dwdx104;
    JSparse[146] = 1.0*dwdx105;
    JSparse[147] = -1.0*dwdx106;
    JSparse[148] = -1.0*dwdx100;
    JSparse[149] = 1.0*dwdx100;
    JSparse[150] = -1.0*dwdx98;
    JSparse[151] = 1.0*dwdx98;
    JSparse[152] = -1.0*dwdx99;
    JSparse[153] = -1.0*dwdx100 - 1.0*dwdx98 - 1.0*dwdx99;
    JSparse[154] = 1.0*dwdx99;
    JSparse[155] = 1.0*dwdx108;
    JSparse[156] = 1.0*dwdx107 + 1.0*dwdx108;
    JSparse[157] = 1.0*dwdx107;
    JSparse[158] = -1.0*dwdx107 - 1.0*dwdx108;
}