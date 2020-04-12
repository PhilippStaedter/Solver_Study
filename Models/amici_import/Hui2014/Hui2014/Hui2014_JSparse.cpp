#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_Hui2014(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0;
    JSparse[1] = 1.0*dwdx1;
    JSparse[2] = -1.0*dwdx2;
    JSparse[3] = 1.0*dwdx3;
    JSparse[4] = -1.0*dwdx3;
    JSparse[5] = 1.0*dwdx3;
    JSparse[6] = 1.0*dwdx4;
    JSparse[7] = -1.0*dwdx5 - 1.0*dwdx6;
    JSparse[8] = 1.0*dwdx5;
    JSparse[9] = -1.0*dwdx5;
    JSparse[10] = 1.0*dwdx7;
    JSparse[11] = -1.0*dwdx7 - 1.0*dwdx8;
    JSparse[12] = 1.0*dwdx7;
    JSparse[13] = -1.0*dwdx8;
    JSparse[14] = 1.0*dwdx8;
    JSparse[15] = -1.0*dwdx10;
    JSparse[16] = 1.0*dwdx10;
    JSparse[17] = -1.0*dwdx10 - 1.0*dwdx11 - 2.0*dwdx9;
    JSparse[18] = 1.0*dwdx9;
    JSparse[19] = 2.0*dwdx12;
    JSparse[20] = -1.0*dwdx12 - 1.0*dwdx13;
    JSparse[21] = -1.0*dwdx13;
    JSparse[22] = 1.0*dwdx13;
    JSparse[23] = -1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17;
    JSparse[24] = 1.0*dwdx15;
    JSparse[25] = 1.0*dwdx16;
    JSparse[26] = 1.0*dwdx17;
    JSparse[27] = -1.0*dwdx15;
    JSparse[28] = -1.0*dwdx16;
    JSparse[29] = -1.0*dwdx17;
    JSparse[30] = 1.0*dwdx14;
    JSparse[31] = -1.0*dwdx14;
    JSparse[32] = 1.0*dwdx18;
    JSparse[33] = -1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    JSparse[34] = 1.0*dwdx19;
    JSparse[35] = 1.0*dwdx20;
    JSparse[36] = 1.0*dwdx18;
    JSparse[37] = -1.0*dwdx19;
    JSparse[38] = -1.0*dwdx20;
    JSparse[39] = 1.0*dwdx21;
    JSparse[40] = 1.0*dwdx22;
    JSparse[41] = -1.0*dwdx21 - 1.0*dwdx22;
    JSparse[42] = 1.0*dwdx21;
    JSparse[43] = 1.0*dwdx22;
    JSparse[44] = 1.0*dwdx23;
    JSparse[45] = 1.0*dwdx24;
    JSparse[46] = -1.0*dwdx23 - 1.0*dwdx24;
    JSparse[47] = 1.0*dwdx23;
    JSparse[48] = 1.0*dwdx24;
    JSparse[49] = -1.0*dwdx29;
    JSparse[50] = 1.0*dwdx29;
    JSparse[51] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    JSparse[52] = 1.0*dwdx30;
    JSparse[53] = 1.0*dwdx31;
    JSparse[54] = -1.0*dwdx30;
    JSparse[55] = -1.0*dwdx31;
    JSparse[56] = -1.0*dwdx25;
    JSparse[57] = 1.0*dwdx25;
    JSparse[58] = -1.0*dwdx34;
    JSparse[59] = 1.0*dwdx34;
    JSparse[60] = 1.0*dwdx33;
    JSparse[61] = -1.0*dwdx33 - 1.0*dwdx34;
    JSparse[62] = 1.0*dwdx33;
    JSparse[63] = -1.0*dwdx32;
    JSparse[64] = 1.0*dwdx32;
    JSparse[65] = -1.0*dwdx36;
    JSparse[66] = 1.0*dwdx36;
    JSparse[67] = 1.0*dwdx35;
    JSparse[68] = -1.0*dwdx35 - 1.0*dwdx36;
    JSparse[69] = 1.0*dwdx35;
    JSparse[70] = -1.0*dwdx41;
    JSparse[71] = 1.0*dwdx41;
    JSparse[72] = -1.0*dwdx38;
    JSparse[73] = 1.0*dwdx38;
    JSparse[74] = -1.0*dwdx38 - 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41;
    JSparse[75] = 1.0*dwdx39 + 1.0*dwdx40;
    JSparse[76] = 1.0*dwdx37;
    JSparse[77] = -1.0*dwdx37;
    JSparse[78] = -1.0*dwdx44;
    JSparse[79] = 1.0*dwdx44;
    JSparse[80] = -1.0*dwdx43;
    JSparse[81] = 1.0*dwdx43;
    JSparse[82] = -1.0*dwdx43 - 1.0*dwdx44;
    JSparse[83] = 1.0*dwdx42;
    JSparse[84] = -1.0*dwdx42;
    JSparse[85] = -1.0*dwdx48;
    JSparse[86] = -1.0*dwdx49;
    JSparse[87] = 1.0*dwdx49;
    JSparse[88] = -1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47;
    JSparse[89] = 1.0*dwdx45 + 1.0*dwdx46 + 1.0*dwdx47;
    JSparse[90] = 1.0*dwdx50 + 1.0*dwdx51 + 1.0*dwdx52;
    JSparse[91] = -1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    JSparse[92] = -1.0*dwdx53;
    JSparse[93] = 1.0*dwdx54;
    JSparse[94] = -1.0*dwdx55;
    JSparse[95] = 1.0*dwdx56;
    JSparse[96] = -1.0*dwdx57;
    JSparse[97] = 1.0*dwdx57;
    JSparse[98] = -1.0*dwdx57;
    JSparse[99] = -1.0*dwdx58 - 1.0*dwdx59;
    JSparse[100] = 1.0*dwdx58 + 1.0*dwdx59;
    JSparse[101] = 1.0*dwdx64;
    JSparse[102] = -1.0*dwdx60;
    JSparse[103] = -1.0*dwdx61;
    JSparse[104] = 1.0*dwdx60;
    JSparse[105] = -1.0*dwdx65;
    JSparse[106] = 1.0*dwdx65;
    JSparse[107] = 1.0*dwdx62;
    JSparse[108] = 1.0*dwdx63;
    JSparse[109] = -1.0*dwdx67;
    JSparse[110] = -1.0*dwdx66 - 1.0*dwdx68;
    JSparse[111] = 1.0*dwdx66 + 1.0*dwdx68;
    JSparse[112] = 1.0*dwdx69;
    JSparse[113] = -1.0*dwdx69;
    JSparse[114] = -1.0*dwdx70;
    JSparse[115] = 1.0*dwdx71;
    JSparse[116] = -1.0*dwdx71;
    JSparse[117] = -1.0*dwdx72;
    JSparse[118] = 1.0*dwdx73;
    JSparse[119] = -1.0*dwdx73;
    JSparse[120] = 1.0*dwdx74;
    JSparse[121] = -1.0*dwdx75;
    JSparse[122] = 1.0*dwdx75;
    JSparse[123] = -1.0*dwdx75 - 1.0*dwdx76;
    JSparse[124] = 1.0*dwdx76;
    JSparse[125] = 1.0*dwdx79;
    JSparse[126] = 1.0*dwdx78;
    JSparse[127] = 1.0*dwdx81;
    JSparse[128] = -1.0*dwdx81;
    JSparse[129] = 1.0*dwdx77;
    JSparse[130] = 1.0*dwdx80;
    JSparse[131] = -1.0*dwdx82 - 1.0*dwdx83;
    JSparse[132] = 1.0*dwdx82 + 1.0*dwdx83;
    JSparse[133] = 1.0*dwdx84;
    JSparse[134] = -1.0*dwdx84;
    JSparse[135] = -1.0*dwdx86;
    JSparse[136] = 1.0*dwdx86;
    JSparse[137] = 1.0*dwdx85;
    JSparse[138] = -1.0*dwdx85;
    JSparse[139] = 1.0*dwdx87;
    JSparse[140] = 1.0*dwdx88;
    JSparse[141] = -1.0*dwdx88;
    JSparse[142] = 1.0*dwdx89;
    JSparse[143] = -1.0*dwdx89;
    JSparse[144] = -1.0*dwdx91;
    JSparse[145] = 1.0*dwdx90;
    JSparse[146] = -1.0*dwdx92;
    JSparse[147] = 1.0*dwdx94;
    JSparse[148] = -1.0*dwdx95;
    JSparse[149] = -1.0*dwdx98;
    JSparse[150] = 1.0*dwdx98;
    JSparse[151] = 1.0*dwdx95;
    JSparse[152] = -1.0*dwdx97;
    JSparse[153] = 1.0*dwdx97;
    JSparse[154] = -1.0*dwdx93 - 1.0*dwdx96;
    JSparse[155] = 1.0*dwdx100;
    JSparse[156] = -1.0*dwdx99;
    JSparse[157] = 1.0*dwdx99;
    JSparse[158] = 1.0*dwdx101;
    JSparse[159] = -1.0*dwdx101;
    JSparse[160] = -1.0*dwdx102;
    JSparse[161] = 1.0*dwdx102;
    JSparse[162] = 1.0*dwdx103 + 1.0*dwdx104;
    JSparse[163] = -1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105;
    JSparse[164] = 1.0*dwdx105;
    JSparse[165] = -1.0*dwdx105;
    JSparse[166] = 1.0*dwdx107;
    JSparse[167] = -1.0*dwdx107;
    JSparse[168] = 1.0*dwdx106;
    JSparse[169] = -1.0*dwdx106;
    JSparse[170] = 1.0*dwdx106;
    JSparse[171] = -1.0*dwdx108;
    JSparse[172] = 1.0*dwdx108;
    JSparse[173] = 1.0*dwdx110;
    JSparse[174] = -1.0*dwdx109 - 1.0*dwdx110;
    JSparse[175] = 1.0*dwdx109;
    JSparse[176] = -1.0*dwdx109;
    JSparse[177] = 1.0*dwdx114;
    JSparse[178] = -1.0*dwdx115;
    JSparse[179] = 1.0*dwdx115;
    JSparse[180] = 1.0*dwdx111;
    JSparse[181] = -1.0*dwdx111;
    JSparse[182] = 1.0*dwdx111;
    JSparse[183] = 1.0*dwdx112;
    JSparse[184] = -1.0*dwdx113;
    JSparse[185] = 1.0*dwdx113;
    JSparse[186] = -1.0*dwdx117;
    JSparse[187] = 1.0*dwdx117;
    JSparse[188] = -1.0*dwdx116;
    JSparse[189] = 1.0*dwdx116;
    JSparse[190] = -1.0*dwdx116 - 1.0*dwdx117;
    JSparse[191] = 1.0*dwdx119;
    JSparse[192] = -1.0*dwdx119;
    JSparse[193] = -1.0*dwdx118 - 1.0*dwdx120 - 1.0*dwdx121;
    JSparse[194] = -1.0*dwdx120;
    JSparse[195] = 1.0*dwdx120;
    JSparse[196] = -1.0*dwdx118;
    JSparse[197] = 1.0*dwdx118;
    JSparse[198] = -1.0*dwdx123;
    JSparse[199] = -1.0*dwdx122;
    JSparse[200] = -1.0*dwdx124 - 1.0*dwdx125;
    JSparse[201] = 1.0*dwdx124;
    JSparse[202] = 1.0*dwdx129;
    JSparse[203] = 1.0*dwdx128;
    JSparse[204] = 1.0*dwdx126;
    JSparse[205] = -1.0*dwdx126;
    JSparse[206] = 1.0*dwdx127;
    JSparse[207] = 1.0*dwdx131;
    JSparse[208] = -1.0*dwdx130;
    JSparse[209] = -1.0*dwdx134;
    JSparse[210] = -1.0*dwdx133;
    JSparse[211] = -1.0*dwdx132 - 1.0*dwdx133 - 1.0*dwdx134;
    JSparse[212] = 1.0*dwdx134;
    JSparse[213] = 1.0*dwdx133;
    JSparse[214] = 1.0*dwdx132;
    JSparse[215] = 1.0*dwdx135;
    JSparse[216] = -1.0*dwdx136;
    JSparse[217] = 1.0*dwdx136;
    JSparse[218] = -1.0*dwdx137;
    JSparse[219] = 1.0*dwdx135;
    JSparse[220] = -1.0*dwdx135 - 1.0*dwdx137;
    JSparse[221] = 1.0*dwdx137;
    JSparse[222] = 1.0*dwdx138;
    JSparse[223] = 1.0*dwdx138;
    JSparse[224] = -1.0*dwdx138 - 1.0*dwdx139;
    JSparse[225] = 1.0*dwdx139;
    JSparse[226] = 1.0*dwdx140;
    JSparse[227] = -1.0*dwdx142;
    JSparse[228] = 1.0*dwdx142;
    JSparse[229] = -1.0*dwdx141;
    JSparse[230] = 1.0*dwdx140;
    JSparse[231] = -1.0*dwdx140 - 1.0*dwdx141;
    JSparse[232] = 1.0*dwdx141;
    JSparse[233] = 1.0*dwdx143;
    JSparse[234] = 1.0*dwdx143;
    JSparse[235] = -1.0*dwdx143 - 1.0*dwdx144;
    JSparse[236] = 1.0*dwdx144;
    JSparse[237] = -1.0*dwdx145;
    JSparse[238] = 1.0*dwdx145;
    JSparse[239] = -1.0*dwdx145;
    JSparse[240] = 1.0*dwdx146;
    JSparse[241] = -1.0*dwdx146;
    JSparse[242] = 1.0*dwdx146;
    JSparse[243] = -1.0*dwdx148;
    JSparse[244] = 1.0*dwdx148;
    JSparse[245] = 1.0*dwdx147;
    JSparse[246] = -1.0*dwdx147 - 1.0*dwdx148;
    JSparse[247] = 1.0*dwdx150;
    JSparse[248] = -1.0*dwdx149;
    JSparse[249] = -1.0*dwdx150;
    JSparse[250] = 1.0*dwdx151 + 1.0*dwdx152;
    JSparse[251] = -1.0*dwdx151 - 1.0*dwdx152;
    JSparse[252] = 1.0*dwdx155;
    JSparse[253] = 1.0*dwdx160;
    JSparse[254] = 1.0*dwdx158;
    JSparse[255] = 1.0*dwdx153;
    JSparse[256] = 1.0*dwdx154;
    JSparse[257] = 1.0*dwdx159;
    JSparse[258] = 1.0*dwdx157;
    JSparse[259] = 1.0*dwdx157;
}