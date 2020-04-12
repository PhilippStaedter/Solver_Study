#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_Liu2011(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0 - 1.0*dwdx1;
    JSparse[1] = -1.0*dwdx0;
    JSparse[2] = 1.0*dwdx0;
    JSparse[3] = -1.0*dwdx1;
    JSparse[4] = 1.0*dwdx1;
    JSparse[5] = -1.0*dwdx2;
    JSparse[6] = -1.0*dwdx2;
    JSparse[7] = 1.0*dwdx2;
    JSparse[8] = -1.0*dwdx3;
    JSparse[9] = -1.0*dwdx3;
    JSparse[10] = 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    JSparse[11] = -1.0*dwdx4;
    JSparse[12] = 1.0*dwdx4;
    JSparse[13] = -1.0*dwdx5;
    JSparse[14] = 1.0*dwdx5;
    JSparse[15] = -1.0*dwdx6;
    JSparse[16] = 1.0*dwdx6;
    JSparse[17] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JSparse[18] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[19] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[20] = -1.0*dwdx15 - 1.0*dwdx16;
    JSparse[21] = -1.0*dwdx15;
    JSparse[22] = 1.0*dwdx15;
    JSparse[23] = -1.0*dwdx16;
    JSparse[24] = 1.0*dwdx16;
    JSparse[25] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JSparse[26] = 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JSparse[27] = 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JSparse[28] = -1.0*dwdx25;
    JSparse[29] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JSparse[30] = 1.0*dwdx25;
    JSparse[31] = -1.0*dwdx27;
    JSparse[32] = -1.0*dwdx28;
    JSparse[33] = -1.0*dwdx26;
    JSparse[34] = 1.0*dwdx26;
    JSparse[35] = 1.0*dwdx27;
    JSparse[36] = 1.0*dwdx28;
    JSparse[37] = -1.0*dwdx29;
    JSparse[38] = -1.0*dwdx30;
    JSparse[39] = 1.0*dwdx30;
    JSparse[40] = 1.0*dwdx30;
    JSparse[41] = -1.0*dwdx31;
    JSparse[42] = -1.0*dwdx29;
    JSparse[43] = 1.0*dwdx29;
    JSparse[44] = 1.0*dwdx31;
    JSparse[45] = 1.0*dwdx31;
    JSparse[46] = -1.0*dwdx32;
    JSparse[47] = -1.0*dwdx32;
    JSparse[48] = 1.0*dwdx32;
    JSparse[49] = -1.0*dwdx33 + 1.0*dwdx37;
    JSparse[50] = -1.0*dwdx33 + 1.0*dwdx37;
    JSparse[51] = 1.0*dwdx33 - 1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38 - 1.0*dwdx39;
    JSparse[52] = -1.0*dwdx34;
    JSparse[53] = 1.0*dwdx34;
    JSparse[54] = 1.0*dwdx34;
    JSparse[55] = 1.0*dwdx35;
    JSparse[56] = -1.0*dwdx38;
    JSparse[57] = 1.0*dwdx36;
    JSparse[58] = 1.0*dwdx38;
    JSparse[59] = -1.0*dwdx40;
    JSparse[60] = 1.0*dwdx40;
    JSparse[61] = 1.0*dwdx40;
    JSparse[62] = -1.0*dwdx41 - 1.0*dwdx42;
    JSparse[63] = 1.0*dwdx41;
    JSparse[64] = -1.0*dwdx43;
    JSparse[65] = 1.0*dwdx43;
    JSparse[66] = -1.0*dwdx44 - 1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47 - 1.0*dwdx48;
    JSparse[67] = -1.0*dwdx45;
    JSparse[68] = 1.0*dwdx45;
    JSparse[69] = -1.0*dwdx44;
    JSparse[70] = 1.0*dwdx44;
    JSparse[71] = -1.0*dwdx46;
    JSparse[72] = -1.0*dwdx47;
    JSparse[73] = 1.0*dwdx46;
    JSparse[74] = 1.0*dwdx47;
    JSparse[75] = -1.0*dwdx48;
    JSparse[76] = 1.0*dwdx48;
    JSparse[77] = -1.0*dwdx49;
    JSparse[78] = 1.0*dwdx49 - 1.0*dwdx50;
    JSparse[79] = -1.0*dwdx50;
    JSparse[80] = 1.0*dwdx50;
    JSparse[81] = -1.0*dwdx51;
    JSparse[82] = 1.0*dwdx51;
    JSparse[83] = -1.0*dwdx51;
    JSparse[84] = -1.0*dwdx54;
    JSparse[85] = -1.0*dwdx53;
    JSparse[86] = -1.0*dwdx52;
    JSparse[87] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JSparse[88] = -1.0*dwdx52;
    JSparse[89] = 1.0*dwdx53;
    JSparse[90] = 1.0*dwdx54;
    JSparse[91] = -1.0*dwdx55;
    JSparse[92] = -1.0*dwdx56;
    JSparse[93] = 1.0*dwdx56;
    JSparse[94] = -1.0*dwdx55 - 1.0*dwdx56;
    JSparse[95] = 1.0*dwdx55;
    JSparse[96] = -1.0*dwdx58;
    JSparse[97] = 1.0*dwdx58;
    JSparse[98] = 1.0*dwdx58;
    JSparse[99] = -1.0*dwdx59;
    JSparse[100] = 1.0*dwdx59;
    JSparse[101] = 1.0*dwdx59;
    JSparse[102] = -1.0*dwdx57;
    JSparse[103] = -1.0*dwdx57;
    JSparse[104] = 1.0*dwdx57;
    JSparse[105] = -1.0*dwdx60;
    JSparse[106] = -1.0*dwdx62;
    JSparse[107] = -1.0*dwdx61;
    JSparse[108] = -1.0*dwdx60;
    JSparse[109] = 1.0*dwdx60 - 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63;
    JSparse[110] = 1.0*dwdx61;
    JSparse[111] = -1.0*dwdx63;
    JSparse[112] = 1.0*dwdx62;
    JSparse[113] = 1.0*dwdx63;
    JSparse[114] = -1.0*dwdx65;
    JSparse[115] = 1.0*dwdx65;
    JSparse[116] = 1.0*dwdx65;
    JSparse[117] = -1.0*dwdx66;
    JSparse[118] = -1.0*dwdx67;
    JSparse[119] = 1.0*dwdx66;
    JSparse[120] = 1.0*dwdx66;
    JSparse[121] = -1.0*dwdx64;
    JSparse[122] = -1.0*dwdx64;
    JSparse[123] = 1.0*dwdx64 - 1.0*dwdx67;
    JSparse[124] = 1.0*dwdx67;
    JSparse[125] = -1.0*dwdx68;
    JSparse[126] = -1.0*dwdx69;
    JSparse[127] = -1.0*dwdx71;
    JSparse[128] = -1.0*dwdx68;
    JSparse[129] = 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71;
    JSparse[130] = 1.0*dwdx69;
    JSparse[131] = -1.0*dwdx70;
    JSparse[132] = 1.0*dwdx70;
    JSparse[133] = 1.0*dwdx71;
    JSparse[134] = -1.0*dwdx73;
    JSparse[135] = 1.0*dwdx73;
    JSparse[136] = 1.0*dwdx73;
    JSparse[137] = -1.0*dwdx74;
    JSparse[138] = -1.0*dwdx72;
    JSparse[139] = 1.0*dwdx74;
    JSparse[140] = 1.0*dwdx74;
    JSparse[141] = -1.0*dwdx72;
    JSparse[142] = 1.0*dwdx72;
    JSparse[143] = -1.0*dwdx75;
    JSparse[144] = -1.0*dwdx78 + 1.0*dwdx79;
    JSparse[145] = 1.0*dwdx79;
    JSparse[146] = -1.0*dwdx77 - 1.0*dwdx79 - 1.0*dwdx80;
    JSparse[147] = -1.0*dwdx81;
    JSparse[148] = -1.0*dwdx83;
    JSparse[149] = -1.0*dwdx76;
    JSparse[150] = -1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx78 - 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83;
    JSparse[151] = 1.0*dwdx75;
    JSparse[152] = 1.0*dwdx76;
    JSparse[153] = 1.0*dwdx77;
    JSparse[154] = 1.0*dwdx78;
    JSparse[155] = 1.0*dwdx80;
    JSparse[156] = 1.0*dwdx81;
    JSparse[157] = 1.0*dwdx83;
    JSparse[158] = -1.0*dwdx84;
    JSparse[159] = -1.0*dwdx84;
    JSparse[160] = 1.0*dwdx84;
    JSparse[161] = -1.0*dwdx85;
    JSparse[162] = -1.0*dwdx85;
    JSparse[163] = 1.0*dwdx85;
    JSparse[164] = -1.0*dwdx86;
    JSparse[165] = -1.0*dwdx86;
    JSparse[166] = 1.0*dwdx86;
    JSparse[167] = -1.0*dwdx87;
    JSparse[168] = -1.0*dwdx87;
    JSparse[169] = 1.0*dwdx87;
    JSparse[170] = -1.0*dwdx88;
    JSparse[171] = -1.0*dwdx88;
    JSparse[172] = 1.0*dwdx88;
    JSparse[173] = -1.0*dwdx90;
    JSparse[174] = 1.0*dwdx90;
    JSparse[175] = 1.0*dwdx90;
    JSparse[176] = -1.0*dwdx91;
    JSparse[177] = -1.0*dwdx89;
    JSparse[178] = 1.0*dwdx91;
    JSparse[179] = 1.0*dwdx91;
    JSparse[180] = -1.0*dwdx92;
    JSparse[181] = -1.0*dwdx89;
    JSparse[182] = 1.0*dwdx89 - 1.0*dwdx92;
    JSparse[183] = 1.0*dwdx92;
    JSparse[184] = -1.0*dwdx93;
    JSparse[185] = -1.0*dwdx93;
    JSparse[186] = 1.0*dwdx93;
    JSparse[187] = -1.0*dwdx95;
    JSparse[188] = 1.0*dwdx95;
    JSparse[189] = 1.0*dwdx95;
    JSparse[190] = -1.0*dwdx96;
    JSparse[191] = 1.0*dwdx96;
    JSparse[192] = 1.0*dwdx96;
    JSparse[193] = -1.0*dwdx94;
    JSparse[194] = -1.0*dwdx94;
    JSparse[195] = 1.0*dwdx94;
    JSparse[196] = -1.0*dwdx99;
    JSparse[197] = 1.0*dwdx99;
    JSparse[198] = 1.0*dwdx99;
    JSparse[199] = -1.0*dwdx100;
    JSparse[200] = -1.0*dwdx98;
    JSparse[201] = 1.0*dwdx100;
    JSparse[202] = 1.0*dwdx100;
    JSparse[203] = -1.0*dwdx97;
    JSparse[204] = -1.0*dwdx98;
    JSparse[205] = -1.0*dwdx97;
    JSparse[206] = 1.0*dwdx97 + 1.0*dwdx98;
    JSparse[207] = -1.0*dwdx102;
    JSparse[208] = 1.0*dwdx101 - 1.0*dwdx102;
    JSparse[209] = -1.0*dwdx101;
    JSparse[210] = 1.0*dwdx102;
    JSparse[211] = -1.0*dwdx101;
    JSparse[212] = 1.0*dwdx103;
    JSparse[213] = -1.0*dwdx103;
    JSparse[214] = -1.0*dwdx103;
    JSparse[215] = -1.0*dwdx105;
    JSparse[216] = 1.0*dwdx105;
    JSparse[217] = 1.0*dwdx105;
    JSparse[218] = -1.0*dwdx106;
    JSparse[219] = 1.0*dwdx106;
    JSparse[220] = 1.0*dwdx106;
    JSparse[221] = -1.0*dwdx104;
    JSparse[222] = -1.0*dwdx104;
    JSparse[223] = 1.0*dwdx104;
    JSparse[224] = 1.0*dwdx107;
    JSparse[225] = -1.0*dwdx107;
    JSparse[226] = -1.0*dwdx107;
}