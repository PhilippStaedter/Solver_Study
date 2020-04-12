#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_jiang1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx1 - 1.0*dwdx2;
    JSparseB[1] = -1.0*dwdx4;
    JSparseB[2] = -1.0*dwdx9;
    JSparseB[3] = -1.0*dwdx23;
    JSparseB[4] = -1.0*dwdx105;
    JSparseB[5] = -1.0*dwdx114;
    JSparseB[6] = -2.0*dwdx3 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6;
    JSparseB[7] = -2.0*dwdx7 - 1.0*dwdx8;
    JSparseB[8] = -2.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 + 1.0*dwdx14;
    JSparseB[9] = 1.0*dwdx35;
    JSparseB[10] = -1.0*dwdx40;
    JSparseB[11] = -1.0*dwdx52;
    JSparseB[12] = 1.0*dwdx109 + 1.0*dwdx110;
    JSparseB[13] = 1.0*dwdx3;
    JSparseB[14] = 1.0*dwdx7;
    JSparseB[15] = 1.0*dwdx11;
    JSparseB[16] = -1.0*dwdx1 + 1.0*dwdx2;
    JSparseB[17] = 1.0*dwdx10 + 1.0*dwdx9;
    JSparseB[18] = 1.0*dwdx23;
    JSparseB[19] = 1.0*dwdx105;
    JSparseB[20] = 1.0*dwdx114;
    JSparseB[21] = 1.0*dwdx3 - 1.0*dwdx5 - 1.0*dwdx6;
    JSparseB[22] = 1.0*dwdx7 + 1.0*dwdx8;
    JSparseB[23] = -1.0*dwdx10;
    JSparseB[24] = 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 - 1.0*dwdx14;
    JSparseB[25] = -1.0*dwdx35;
    JSparseB[26] = 1.0*dwdx40;
    JSparseB[27] = 1.0*dwdx52;
    JSparseB[28] = -1.0*dwdx109 - 1.0*dwdx110;
    JSparseB[29] = 1.0*dwdx15 - 1.0*dwdx16;
    JSparseB[30] = -1.0*dwdx31;
    JSparseB[31] = -1.0*dwdx84;
    JSparseB[32] = -1.0*dwdx94;
    JSparseB[33] = 1.0*dwdx102;
    JSparseB[34] = -1.0*dwdx115;
    JSparseB[35] = -1.0*dwdx17;
    JSparseB[36] = -1.0*dwdx106;
    JSparseB[37] = 1.0*dwdx18;
    JSparseB[38] = 1.0*dwdx54;
    JSparseB[39] = 1.0*dwdx96;
    JSparseB[40] = 1.0*dwdx113;
    JSparseB[41] = -1.0*dwdx19 + 1.0*dwdx20;
    JSparseB[42] = 1.0*dwdx21;
    JSparseB[43] = -1.0*dwdx55 + 1.0*dwdx56;
    JSparseB[44] = 1.0*dwdx57;
    JSparseB[45] = -1.0*dwdx97;
    JSparseB[46] = -1.0*dwdx104;
    JSparseB[47] = -1.0*dwdx20;
    JSparseB[48] = -1.0*dwdx21 + 1.0*dwdx22;
    JSparseB[49] = -1.0*dwdx56;
    JSparseB[50] = -1.0*dwdx57 + 1.0*dwdx58;
    JSparseB[51] = 1.0*dwdx100;
    JSparseB[52] = 1.0*dwdx108;
    JSparseB[53] = 1.0*dwdx2;
    JSparseB[54] = 1.0*dwdx9;
    JSparseB[55] = -1.0*dwdx16;
    JSparseB[56] = 1.0*dwdx23;
    JSparseB[57] = -1.0*dwdx29 - 1.0*dwdx31;
    JSparseB[58] = -1.0*dwdx81 - 1.0*dwdx84;
    JSparseB[59] = -1.0*dwdx91 - 1.0*dwdx94;
    JSparseB[60] = -1.0*dwdx95;
    JSparseB[61] = 1.0*dwdx105;
    JSparseB[62] = 1.0*dwdx114 - 1.0*dwdx115;
    JSparseB[63] = -1.0*dwdx123;
    JSparseB[64] = -1.0*dwdx15;
    JSparseB[65] = 1.0*dwdx25 - 1.0*dwdx26;
    JSparseB[66] = -1.0*dwdx28;
    JSparseB[67] = 1.0*dwdx59;
    JSparseB[68] = -1.0*dwdx70;
    JSparseB[69] = -1.0*dwdx75;
    JSparseB[70] = -1.0*dwdx102;
    JSparseB[71] = 1.0*dwdx17;
    JSparseB[72] = 1.0*dwdx26;
    JSparseB[73] = 1.0*dwdx27 + 1.0*dwdx28;
    JSparseB[74] = 1.0*dwdx62;
    JSparseB[75] = 1.0*dwdx70;
    JSparseB[76] = 1.0*dwdx75;
    JSparseB[77] = 1.0*dwdx106;
    JSparseB[78] = -1.0*dwdx15 + 1.0*dwdx16;
    JSparseB[79] = 1.0*dwdx29 - 1.0*dwdx30 + 1.0*dwdx31;
    JSparseB[80] = -1.0*dwdx50;
    JSparseB[81] = -1.0*dwdx53;
    JSparseB[82] = 1.0*dwdx81 + 1.0*dwdx84;
    JSparseB[83] = 1.0*dwdx91 + 1.0*dwdx94;
    JSparseB[84] = 1.0*dwdx95;
    JSparseB[85] = -1.0*dwdx102;
    JSparseB[86] = 1.0*dwdx115;
    JSparseB[87] = 1.0*dwdx123 - 1.0*dwdx124;
    JSparseB[88] = -1.0*dwdx125;
    JSparseB[89] = 1.0*dwdx17;
    JSparseB[90] = 1.0*dwdx106;
    JSparseB[91] = -2.0*dwdx32 + 1.0*dwdx33;
    JSparseB[92] = -2.0*dwdx34;
    JSparseB[93] = -2.0*dwdx121;
    JSparseB[94] = 2.0*dwdx32 - 1.0*dwdx33;
    JSparseB[95] = 2.0*dwdx34;
    JSparseB[96] = 2.0*dwdx121;
    JSparseB[97] = -1.0*dwdx46;
    JSparseB[98] = 1.0*dwdx86;
    JSparseB[99] = 1.0*dwdx5;
    JSparseB[100] = 1.0*dwdx14;
    JSparseB[101] = 1.0*dwdx35;
    JSparseB[102] = -1.0*dwdx49;
    JSparseB[103] = -1.0*dwdx79;
    JSparseB[104] = 1.0*dwdx109;
    JSparseB[105] = -1.0*dwdx36 + 1.0*dwdx37;
    JSparseB[106] = -1.0*dwdx38 + 1.0*dwdx39;
    JSparseB[107] = 1.0*dwdx41;
    JSparseB[108] = 1.0*dwdx42;
    JSparseB[109] = -1.0*dwdx118;
    JSparseB[110] = -1.0*dwdx122;
    JSparseB[111] = 1.0*dwdx36 - 1.0*dwdx37;
    JSparseB[112] = 1.0*dwdx38 - 1.0*dwdx39;
    JSparseB[113] = -1.0*dwdx41;
    JSparseB[114] = -1.0*dwdx42;
    JSparseB[115] = 1.0*dwdx118;
    JSparseB[116] = 1.0*dwdx122;
    JSparseB[117] = 1.0*dwdx8;
    JSparseB[118] = -1.0*dwdx12 + 1.0*dwdx13;
    JSparseB[119] = 1.0*dwdx40;
    JSparseB[120] = -1.0*dwdx52;
    JSparseB[121] = -1.0*dwdx37;
    JSparseB[122] = -1.0*dwdx39;
    JSparseB[123] = -1.0*dwdx41;
    JSparseB[124] = -1.0*dwdx42;
    JSparseB[125] = 1.0*dwdx46;
    JSparseB[126] = 1.0*dwdx37;
    JSparseB[127] = 1.0*dwdx39;
    JSparseB[128] = 1.0*dwdx41;
    JSparseB[129] = 1.0*dwdx42;
    JSparseB[130] = -1.0*dwdx46;
    JSparseB[131] = -1.0*dwdx8;
    JSparseB[132] = -1.0*dwdx13;
    JSparseB[133] = -1.0*dwdx40;
    JSparseB[134] = 1.0*dwdx43;
    JSparseB[135] = 1.0*dwdx48;
    JSparseB[136] = -1.0*dwdx44 + 1.0*dwdx45;
    JSparseB[137] = 1.0*dwdx67;
    JSparseB[138] = -1.0*dwdx116;
    JSparseB[139] = -1.0*dwdx119;
    JSparseB[140] = -1.0*dwdx126;
    JSparseB[141] = 1.0*dwdx46;
    JSparseB[142] = -1.0*dwdx86;
    JSparseB[143] = -2.0*dwdx43;
    JSparseB[144] = 1.0*dwdx47 - 2.0*dwdx48 + 1.0*dwdx49;
    JSparseB[145] = 1.0*dwdx79;
    JSparseB[146] = 1.0*dwdx30;
    JSparseB[147] = 1.0*dwdx50;
    JSparseB[148] = 1.0*dwdx53;
    JSparseB[149] = 1.0*dwdx124;
    JSparseB[150] = 1.0*dwdx125;
    JSparseB[151] = 1.0*dwdx12;
    JSparseB[152] = -1.0*dwdx51 + 1.0*dwdx52;
    JSparseB[153] = -1.0*dwdx30;
    JSparseB[154] = -1.0*dwdx50;
    JSparseB[155] = -1.0*dwdx53;
    JSparseB[156] = -1.0*dwdx124;
    JSparseB[157] = -1.0*dwdx125;
    JSparseB[158] = -1.0*dwdx18;
    JSparseB[159] = 1.0*dwdx19 - 1.0*dwdx20;
    JSparseB[160] = -1.0*dwdx21;
    JSparseB[161] = -1.0*dwdx54 + 1.0*dwdx55 - 1.0*dwdx56;
    JSparseB[162] = -1.0*dwdx57;
    JSparseB[163] = -1.0*dwdx96 + 1.0*dwdx97;
    JSparseB[164] = 1.0*dwdx104;
    JSparseB[165] = -1.0*dwdx113;
    JSparseB[166] = 1.0*dwdx20;
    JSparseB[167] = 1.0*dwdx21 - 1.0*dwdx22;
    JSparseB[168] = 1.0*dwdx56;
    JSparseB[169] = 1.0*dwdx57 - 1.0*dwdx58;
    JSparseB[170] = -1.0*dwdx100;
    JSparseB[171] = -1.0*dwdx108;
    JSparseB[172] = -1.0*dwdx1;
    JSparseB[173] = 1.0*dwdx0;
    JSparseB[174] = -1.0*dwdx25;
    JSparseB[175] = -1.0*dwdx59 + 1.0*dwdx60 - 1.0*dwdx61;
    JSparseB[176] = -1.0*dwdx64;
    JSparseB[177] = -1.0*dwdx71;
    JSparseB[178] = -1.0*dwdx77;
    JSparseB[179] = 1.0*dwdx24;
    JSparseB[180] = -1.0*dwdx27;
    JSparseB[181] = 1.0*dwdx61;
    JSparseB[182] = -1.0*dwdx62 + 1.0*dwdx63 + 1.0*dwdx64;
    JSparseB[183] = 1.0*dwdx71;
    JSparseB[184] = 1.0*dwdx77;
    JSparseB[185] = 1.0*dwdx88;
    JSparseB[186] = 1.0*dwdx90;
    JSparseB[187] = 1.0*dwdx101;
    JSparseB[188] = 1.0*dwdx65 - 1.0*dwdx66;
    JSparseB[189] = -1.0*dwdx80;
    JSparseB[190] = -1.0*dwdx87;
    JSparseB[191] = -1.0*dwdx111;
    JSparseB[192] = 1.0*dwdx26;
    JSparseB[193] = 1.0*dwdx28;
    JSparseB[194] = -1.0*dwdx45;
    JSparseB[195] = 1.0*dwdx61;
    JSparseB[196] = 1.0*dwdx64;
    JSparseB[197] = -1.0*dwdx67 + 1.0*dwdx68 - 1.0*dwdx69 + 1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72;
    JSparseB[198] = -1.0*dwdx73 + 1.0*dwdx75 + 1.0*dwdx77;
    JSparseB[199] = 1.0*dwdx82;
    JSparseB[200] = 1.0*dwdx92;
    JSparseB[201] = -1.0*dwdx98;
    JSparseB[202] = -1.0*dwdx99;
    JSparseB[203] = 1.0*dwdx103;
    JSparseB[204] = -1.0*dwdx26;
    JSparseB[205] = -1.0*dwdx28;
    JSparseB[206] = -1.0*dwdx61;
    JSparseB[207] = -1.0*dwdx64;
    JSparseB[208] = 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71;
    JSparseB[209] = 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75 + 1.0*dwdx76 - 1.0*dwdx77;
    JSparseB[210] = -1.0*dwdx78;
    JSparseB[211] = -1.0*dwdx85;
    JSparseB[212] = 1.0*dwdx89;
    JSparseB[213] = 1.0*dwdx98;
    JSparseB[214] = 1.0*dwdx99;
    JSparseB[215] = -1.0*dwdx107;
    JSparseB[216] = 1.0*dwdx49;
    JSparseB[217] = -1.0*dwdx66;
    JSparseB[218] = -1.0*dwdx74;
    JSparseB[219] = -1.0*dwdx78 + 1.0*dwdx79 - 1.0*dwdx80;
    JSparseB[220] = -1.0*dwdx85 - 1.0*dwdx86 - 1.0*dwdx87;
    JSparseB[221] = -1.0*dwdx107;
    JSparseB[222] = -1.0*dwdx111;
    JSparseB[223] = -1.0*dwdx0;
    JSparseB[224] = -1.0*dwdx16;
    JSparseB[225] = -1.0*dwdx29 - 1.0*dwdx31;
    JSparseB[226] = -1.0*dwdx60;
    JSparseB[227] = -1.0*dwdx68;
    JSparseB[228] = -1.0*dwdx81 - 1.0*dwdx82 + 1.0*dwdx83 - 1.0*dwdx84;
    JSparseB[229] = -1.0*dwdx91 - 1.0*dwdx92 + 1.0*dwdx93 - 1.0*dwdx94;
    JSparseB[230] = -1.0*dwdx95;
    JSparseB[231] = -1.0*dwdx103;
    JSparseB[232] = -1.0*dwdx115;
    JSparseB[233] = 1.0*dwdx117;
    JSparseB[234] = 1.0*dwdx120;
    JSparseB[235] = -1.0*dwdx123;
    JSparseB[236] = -1.0*dwdx49;
    JSparseB[237] = 1.0*dwdx66;
    JSparseB[238] = 1.0*dwdx74;
    JSparseB[239] = 1.0*dwdx78 - 1.0*dwdx79 + 1.0*dwdx80;
    JSparseB[240] = 1.0*dwdx85 + 1.0*dwdx86 + 1.0*dwdx87;
    JSparseB[241] = 1.0*dwdx107;
    JSparseB[242] = 1.0*dwdx111;
    JSparseB[243] = -1.0*dwdx72;
    JSparseB[244] = -1.0*dwdx24;
    JSparseB[245] = -1.0*dwdx63;
    JSparseB[246] = -1.0*dwdx76;
    JSparseB[247] = -1.0*dwdx88;
    JSparseB[248] = -1.0*dwdx89 - 1.0*dwdx90;
    JSparseB[249] = -1.0*dwdx101;
    JSparseB[250] = 1.0*dwdx24;
    JSparseB[251] = 1.0*dwdx63;
    JSparseB[252] = 1.0*dwdx76;
    JSparseB[253] = 1.0*dwdx88;
    JSparseB[254] = 1.0*dwdx89 + 1.0*dwdx90;
    JSparseB[255] = 1.0*dwdx101;
    JSparseB[256] = 1.0*dwdx72;
    JSparseB[257] = 1.0*dwdx0;
    JSparseB[258] = 1.0*dwdx16;
    JSparseB[259] = 1.0*dwdx29 + 1.0*dwdx31;
    JSparseB[260] = 1.0*dwdx60;
    JSparseB[261] = 1.0*dwdx68;
    JSparseB[262] = 1.0*dwdx81 + 1.0*dwdx82 - 1.0*dwdx83 + 1.0*dwdx84;
    JSparseB[263] = 1.0*dwdx91 + 1.0*dwdx92 - 1.0*dwdx93 + 1.0*dwdx94;
    JSparseB[264] = 1.0*dwdx95;
    JSparseB[265] = 1.0*dwdx103;
    JSparseB[266] = 1.0*dwdx115;
    JSparseB[267] = -1.0*dwdx117;
    JSparseB[268] = -1.0*dwdx120;
    JSparseB[269] = 1.0*dwdx123;
    JSparseB[270] = -1.0*dwdx0;
    JSparseB[271] = 1.0*dwdx18;
    JSparseB[272] = -1.0*dwdx19;
    JSparseB[273] = 1.0*dwdx29;
    JSparseB[274] = 1.0*dwdx54 - 1.0*dwdx55;
    JSparseB[275] = -1.0*dwdx60;
    JSparseB[276] = 1.0*dwdx69;
    JSparseB[277] = 1.0*dwdx73;
    JSparseB[278] = 1.0*dwdx81;
    JSparseB[279] = 1.0*dwdx91;
    JSparseB[280] = 1.0*dwdx95 + 1.0*dwdx96 - 1.0*dwdx97 + 1.0*dwdx98;
    JSparseB[281] = 1.0*dwdx99;
    JSparseB[282] = -1.0*dwdx104;
    JSparseB[283] = 1.0*dwdx113;
    JSparseB[284] = 1.0*dwdx123;
    JSparseB[285] = 1.0*dwdx22;
    JSparseB[286] = -1.0*dwdx24;
    JSparseB[287] = 1.0*dwdx58;
    JSparseB[288] = -1.0*dwdx63;
    JSparseB[289] = -1.0*dwdx69;
    JSparseB[290] = -1.0*dwdx73;
    JSparseB[291] = -1.0*dwdx88;
    JSparseB[292] = -1.0*dwdx90;
    JSparseB[293] = -1.0*dwdx98;
    JSparseB[294] = 1.0*dwdx100 - 1.0*dwdx101 - 1.0*dwdx99;
    JSparseB[295] = 1.0*dwdx108;
    JSparseB[296] = -1.0*dwdx2;
    JSparseB[297] = -1.0*dwdx9;
    JSparseB[298] = 1.0*dwdx15;
    JSparseB[299] = 1.0*dwdx19;
    JSparseB[300] = -1.0*dwdx23;
    JSparseB[301] = 1.0*dwdx55;
    JSparseB[302] = -1.0*dwdx68;
    JSparseB[303] = -1.0*dwdx82;
    JSparseB[304] = -1.0*dwdx92;
    JSparseB[305] = 1.0*dwdx97;
    JSparseB[306] = 1.0*dwdx102 - 1.0*dwdx103 + 1.0*dwdx104 - 1.0*dwdx105;
    JSparseB[307] = -1.0*dwdx114;
    JSparseB[308] = -1.0*dwdx17;
    JSparseB[309] = -1.0*dwdx22;
    JSparseB[310] = -1.0*dwdx58;
    JSparseB[311] = 1.0*dwdx74;
    JSparseB[312] = 1.0*dwdx78;
    JSparseB[313] = 1.0*dwdx85;
    JSparseB[314] = -1.0*dwdx100;
    JSparseB[315] = -1.0*dwdx106 + 1.0*dwdx107 - 1.0*dwdx108;
    JSparseB[316] = -1.0*dwdx5 + 1.0*dwdx6;
    JSparseB[317] = -1.0*dwdx14;
    JSparseB[318] = -1.0*dwdx35;
    JSparseB[319] = -1.0*dwdx109 + 1.0*dwdx110;
    JSparseB[320] = -1.0*dwdx6;
    JSparseB[321] = 1.0*dwdx66;
    JSparseB[322] = -1.0*dwdx76;
    JSparseB[323] = 1.0*dwdx80;
    JSparseB[324] = 1.0*dwdx87;
    JSparseB[325] = -1.0*dwdx89;
    JSparseB[326] = -1.0*dwdx110;
    JSparseB[327] = 1.0*dwdx111 + 1.0*dwdx112;
    JSparseB[328] = 1.0*dwdx1 - 1.0*dwdx2;
    JSparseB[329] = -1.0*dwdx9;
    JSparseB[330] = -1.0*dwdx23;
    JSparseB[331] = 1.0*dwdx30;
    JSparseB[332] = 1.0*dwdx50;
    JSparseB[333] = 1.0*dwdx53;
    JSparseB[334] = -1.0*dwdx105;
    JSparseB[335] = -1.0*dwdx114;
    JSparseB[336] = 1.0*dwdx124;
    JSparseB[337] = 1.0*dwdx125;
    JSparseB[338] = 1.0*dwdx2;
    JSparseB[339] = 1.0*dwdx9;
    JSparseB[340] = 1.0*dwdx16;
    JSparseB[341] = -1.0*dwdx18;
    JSparseB[342] = 1.0*dwdx23;
    JSparseB[343] = 1.0*dwdx31;
    JSparseB[344] = -1.0*dwdx54;
    JSparseB[345] = -1.0*dwdx72;
    JSparseB[346] = 1.0*dwdx84;
    JSparseB[347] = 1.0*dwdx94;
    JSparseB[348] = -1.0*dwdx96;
    JSparseB[349] = 1.0*dwdx105;
    JSparseB[350] = -1.0*dwdx112;
    JSparseB[351] = -1.0*dwdx113 + 1.0*dwdx114 + 1.0*dwdx115;
    JSparseB[352] = -1.0*dwdx32;
    JSparseB[353] = -1.0*dwdx34;
    JSparseB[354] = 1.0*dwdx36;
    JSparseB[355] = 1.0*dwdx38;
    JSparseB[356] = 1.0*dwdx44;
    JSparseB[357] = 1.0*dwdx83;
    JSparseB[358] = 1.0*dwdx93;
    JSparseB[359] = 1.0*dwdx116 + 1.0*dwdx117 + 1.0*dwdx118;
    JSparseB[360] = 1.0*dwdx119 + 1.0*dwdx120 - 1.0*dwdx121 + 1.0*dwdx122;
    JSparseB[361] = 1.0*dwdx126;
    JSparseB[362] = 1.0*dwdx32;
    JSparseB[363] = 1.0*dwdx34;
    JSparseB[364] = -1.0*dwdx36;
    JSparseB[365] = -1.0*dwdx38;
    JSparseB[366] = -1.0*dwdx44;
    JSparseB[367] = -1.0*dwdx83;
    JSparseB[368] = -1.0*dwdx93;
    JSparseB[369] = -1.0*dwdx116 - 1.0*dwdx117 - 1.0*dwdx118;
    JSparseB[370] = -1.0*dwdx119 - 1.0*dwdx120 + 1.0*dwdx121 - 1.0*dwdx122;
    JSparseB[371] = -1.0*dwdx126;
    JSparseB[372] = -1.0*dwdx29 + 1.0*dwdx30;
    JSparseB[373] = 1.0*dwdx50;
    JSparseB[374] = 1.0*dwdx53;
    JSparseB[375] = -1.0*dwdx81;
    JSparseB[376] = -1.0*dwdx91;
    JSparseB[377] = -1.0*dwdx95;
    JSparseB[378] = -1.0*dwdx123 + 1.0*dwdx124;
    JSparseB[379] = 1.0*dwdx125;
    JSparseB[380] = -1.0*dwdx30;
    JSparseB[381] = 1.0*dwdx44;
    JSparseB[382] = -1.0*dwdx50;
    JSparseB[383] = -1.0*dwdx53;
    JSparseB[384] = 1.0*dwdx116;
    JSparseB[385] = 1.0*dwdx119;
    JSparseB[386] = -1.0*dwdx124;
    JSparseB[387] = -1.0*dwdx125 + 1.0*dwdx126;
}