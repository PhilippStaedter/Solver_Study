#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_Liu2011(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JSparseB[1] = 1.0*dwdx2;
    JSparseB[2] = 1.0*dwdx3;
    JSparseB[3] = 1.0*dwdx54;
    JSparseB[4] = 1.0*dwdx68;
    JSparseB[5] = 1.0*dwdx0;
    JSparseB[6] = 1.0*dwdx2;
    JSparseB[7] = 1.0*dwdx3;
    JSparseB[8] = -1.0*dwdx0;
    JSparseB[9] = -1.0*dwdx2;
    JSparseB[10] = -1.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JSparseB[11] = 1.0*dwdx25;
    JSparseB[12] = 1.0*dwdx29;
    JSparseB[13] = 1.0*dwdx55;
    JSparseB[14] = 1.0*dwdx60;
    JSparseB[15] = 1.0*dwdx75;
    JSparseB[16] = 1.0*dwdx84;
    JSparseB[17] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JSparseB[18] = 1.0*dwdx30;
    JSparseB[19] = 1.0*dwdx58;
    JSparseB[20] = 1.0*dwdx65;
    JSparseB[21] = 1.0*dwdx73;
    JSparseB[22] = 1.0*dwdx90;
    JSparseB[23] = 1.0*dwdx95;
    JSparseB[24] = 1.0*dwdx99;
    JSparseB[25] = 1.0*dwdx105;
    JSparseB[26] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JSparseB[27] = -1.0*dwdx30;
    JSparseB[28] = -1.0*dwdx58;
    JSparseB[29] = -1.0*dwdx65;
    JSparseB[30] = -1.0*dwdx73;
    JSparseB[31] = -1.0*dwdx90;
    JSparseB[32] = -1.0*dwdx95;
    JSparseB[33] = -1.0*dwdx99;
    JSparseB[34] = -1.0*dwdx105;
    JSparseB[35] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JSparseB[36] = 1.0*dwdx15 + 1.0*dwdx16;
    JSparseB[37] = -1.0*dwdx30;
    JSparseB[38] = 1.0*dwdx32;
    JSparseB[39] = 1.0*dwdx33 - 1.0*dwdx37;
    JSparseB[40] = -1.0*dwdx58;
    JSparseB[41] = -1.0*dwdx65;
    JSparseB[42] = -1.0*dwdx73;
    JSparseB[43] = 1.0*dwdx78 - 1.0*dwdx79;
    JSparseB[44] = 1.0*dwdx86;
    JSparseB[45] = -1.0*dwdx90;
    JSparseB[46] = -1.0*dwdx95;
    JSparseB[47] = -1.0*dwdx99;
    JSparseB[48] = -1.0*dwdx105;
    JSparseB[49] = 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JSparseB[50] = 1.0*dwdx31;
    JSparseB[51] = 1.0*dwdx59;
    JSparseB[52] = 1.0*dwdx66;
    JSparseB[53] = 1.0*dwdx74;
    JSparseB[54] = 1.0*dwdx91;
    JSparseB[55] = 1.0*dwdx96;
    JSparseB[56] = 1.0*dwdx100;
    JSparseB[57] = 1.0*dwdx106;
    JSparseB[58] = 1.0*dwdx4;
    JSparseB[59] = 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28;
    JSparseB[60] = 1.0*dwdx29;
    JSparseB[61] = 1.0*dwdx62;
    JSparseB[62] = 1.0*dwdx67;
    JSparseB[63] = 1.0*dwdx69;
    JSparseB[64] = 1.0*dwdx72;
    JSparseB[65] = 1.0*dwdx89;
    JSparseB[66] = 1.0*dwdx98;
    JSparseB[67] = -1.0*dwdx4;
    JSparseB[68] = -1.0*dwdx25;
    JSparseB[69] = -1.0*dwdx29;
    JSparseB[70] = 1.0*dwdx15;
    JSparseB[71] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JSparseB[72] = -1.0*dwdx31;
    JSparseB[73] = 1.0*dwdx32;
    JSparseB[74] = 1.0*dwdx33 - 1.0*dwdx37;
    JSparseB[75] = -1.0*dwdx59;
    JSparseB[76] = -1.0*dwdx66;
    JSparseB[77] = -1.0*dwdx74;
    JSparseB[78] = -1.0*dwdx79;
    JSparseB[79] = -1.0*dwdx91;
    JSparseB[80] = -1.0*dwdx96;
    JSparseB[81] = -1.0*dwdx100;
    JSparseB[82] = -1.0*dwdx106;
    JSparseB[83] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JSparseB[84] = -1.0*dwdx31;
    JSparseB[85] = -1.0*dwdx59;
    JSparseB[86] = -1.0*dwdx66;
    JSparseB[87] = -1.0*dwdx74;
    JSparseB[88] = -1.0*dwdx91;
    JSparseB[89] = -1.0*dwdx96;
    JSparseB[90] = -1.0*dwdx100;
    JSparseB[91] = -1.0*dwdx106;
    JSparseB[92] = -1.0*dwdx15;
    JSparseB[93] = -1.0*dwdx32;
    JSparseB[94] = -1.0*dwdx33 + 1.0*dwdx35 + 1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx38 + 1.0*dwdx39;
    JSparseB[95] = 1.0*dwdx49;
    JSparseB[96] = 1.0*dwdx77 + 1.0*dwdx79 + 1.0*dwdx80;
    JSparseB[97] = 1.0*dwdx87;
    JSparseB[98] = 1.0*dwdx34;
    JSparseB[99] = 1.0*dwdx40;
    JSparseB[100] = -1.0*dwdx34;
    JSparseB[101] = -1.0*dwdx40;
    JSparseB[102] = -1.0*dwdx34;
    JSparseB[103] = -1.0*dwdx40;
    JSparseB[104] = 1.0*dwdx41 + 1.0*dwdx42;
    JSparseB[105] = 1.0*dwdx43;
    JSparseB[106] = -1.0*dwdx41;
    JSparseB[107] = -1.0*dwdx43;
    JSparseB[108] = 1.0*dwdx44 + 1.0*dwdx45 + 1.0*dwdx46 + 1.0*dwdx47 + 1.0*dwdx48;
    JSparseB[109] = 1.0*dwdx53;
    JSparseB[110] = 1.0*dwdx57;
    JSparseB[111] = 1.0*dwdx61;
    JSparseB[112] = 1.0*dwdx64;
    JSparseB[113] = 1.0*dwdx71;
    JSparseB[114] = 1.0*dwdx92;
    JSparseB[115] = 1.0*dwdx94;
    JSparseB[116] = 1.0*dwdx97;
    JSparseB[117] = 1.0*dwdx102;
    JSparseB[118] = 1.0*dwdx104;
    JSparseB[119] = -1.0*dwdx35;
    JSparseB[120] = -1.0*dwdx49 + 1.0*dwdx50;
    JSparseB[121] = 1.0*dwdx81;
    JSparseB[122] = 1.0*dwdx88;
    JSparseB[123] = 1.0*dwdx51;
    JSparseB[124] = 1.0*dwdx52;
    JSparseB[125] = 1.0*dwdx56;
    JSparseB[126] = 1.0*dwdx1;
    JSparseB[127] = 1.0*dwdx45;
    JSparseB[128] = -1.0*dwdx51;
    JSparseB[129] = -1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54;
    JSparseB[130] = -1.0*dwdx56;
    JSparseB[131] = 1.0*dwdx57;
    JSparseB[132] = 1.0*dwdx68;
    JSparseB[133] = 1.0*dwdx5;
    JSparseB[134] = 1.0*dwdx51;
    JSparseB[135] = 1.0*dwdx52;
    JSparseB[136] = 1.0*dwdx55 + 1.0*dwdx56;
    JSparseB[137] = 1.0*dwdx60;
    JSparseB[138] = -1.0*dwdx45;
    JSparseB[139] = -1.0*dwdx53;
    JSparseB[140] = -1.0*dwdx57;
    JSparseB[141] = -1.0*dwdx5;
    JSparseB[142] = 1.0*dwdx27;
    JSparseB[143] = 1.0*dwdx44;
    JSparseB[144] = -1.0*dwdx55;
    JSparseB[145] = -1.0*dwdx60 + 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63;
    JSparseB[146] = 1.0*dwdx64;
    JSparseB[147] = 1.0*dwdx83;
    JSparseB[148] = 1.0*dwdx89;
    JSparseB[149] = 1.0*dwdx93;
    JSparseB[150] = 1.0*dwdx28;
    JSparseB[151] = -1.0*dwdx44;
    JSparseB[152] = -1.0*dwdx61;
    JSparseB[153] = -1.0*dwdx64 + 1.0*dwdx67;
    JSparseB[154] = 1.0*dwdx98;
    JSparseB[155] = -1.0*dwdx1;
    JSparseB[156] = 1.0*dwdx26;
    JSparseB[157] = 1.0*dwdx46;
    JSparseB[158] = -1.0*dwdx54;
    JSparseB[159] = -1.0*dwdx68 + 1.0*dwdx69 + 1.0*dwdx70 + 1.0*dwdx71;
    JSparseB[160] = 1.0*dwdx72;
    JSparseB[161] = 1.0*dwdx76;
    JSparseB[162] = 1.0*dwdx85;
    JSparseB[163] = 1.0*dwdx94;
    JSparseB[164] = -1.0*dwdx26;
    JSparseB[165] = -1.0*dwdx69;
    JSparseB[166] = -1.0*dwdx72;
    JSparseB[167] = 1.0*dwdx6;
    JSparseB[168] = 1.0*dwdx16;
    JSparseB[169] = 1.0*dwdx38;
    JSparseB[170] = 1.0*dwdx50;
    JSparseB[171] = 1.0*dwdx63;
    JSparseB[172] = 1.0*dwdx70;
    JSparseB[173] = 1.0*dwdx75 + 1.0*dwdx76 + 1.0*dwdx78 + 1.0*dwdx80 + 1.0*dwdx81 + 1.0*dwdx82 + 1.0*dwdx83;
    JSparseB[174] = 1.0*dwdx84;
    JSparseB[175] = 1.0*dwdx85;
    JSparseB[176] = 1.0*dwdx86;
    JSparseB[177] = 1.0*dwdx87;
    JSparseB[178] = 1.0*dwdx88;
    JSparseB[179] = 1.0*dwdx93;
    JSparseB[180] = -1.0*dwdx6;
    JSparseB[181] = -1.0*dwdx75;
    JSparseB[182] = -1.0*dwdx84;
    JSparseB[183] = -1.0*dwdx70;
    JSparseB[184] = -1.0*dwdx76;
    JSparseB[185] = -1.0*dwdx85;
    JSparseB[186] = -1.0*dwdx36;
    JSparseB[187] = -1.0*dwdx77;
    JSparseB[188] = -1.0*dwdx16;
    JSparseB[189] = -1.0*dwdx78;
    JSparseB[190] = -1.0*dwdx86;
    JSparseB[191] = -1.0*dwdx38;
    JSparseB[192] = -1.0*dwdx80;
    JSparseB[193] = -1.0*dwdx87;
    JSparseB[194] = -1.0*dwdx50;
    JSparseB[195] = -1.0*dwdx81;
    JSparseB[196] = -1.0*dwdx88;
    JSparseB[197] = -1.0*dwdx27;
    JSparseB[198] = 1.0*dwdx47;
    JSparseB[199] = -1.0*dwdx62;
    JSparseB[200] = -1.0*dwdx89 + 1.0*dwdx92;
    JSparseB[201] = 1.0*dwdx97;
    JSparseB[202] = -1.0*dwdx63;
    JSparseB[203] = -1.0*dwdx83;
    JSparseB[204] = -1.0*dwdx93;
    JSparseB[205] = -1.0*dwdx46;
    JSparseB[206] = -1.0*dwdx71;
    JSparseB[207] = -1.0*dwdx94;
    JSparseB[208] = -1.0*dwdx28;
    JSparseB[209] = -1.0*dwdx47;
    JSparseB[210] = -1.0*dwdx67;
    JSparseB[211] = -1.0*dwdx92;
    JSparseB[212] = -1.0*dwdx97 - 1.0*dwdx98;
    JSparseB[213] = 1.0*dwdx48;
    JSparseB[214] = -1.0*dwdx101 + 1.0*dwdx102;
    JSparseB[215] = -1.0*dwdx103;
    JSparseB[216] = 1.0*dwdx104;
    JSparseB[217] = -1.0*dwdx107;
    JSparseB[218] = 1.0*dwdx101;
    JSparseB[219] = 1.0*dwdx103;
    JSparseB[220] = 1.0*dwdx107;
    JSparseB[221] = -1.0*dwdx48;
    JSparseB[222] = -1.0*dwdx102;
    JSparseB[223] = -1.0*dwdx104;
    JSparseB[224] = 1.0*dwdx101;
    JSparseB[225] = 1.0*dwdx103;
    JSparseB[226] = 1.0*dwdx107;
}