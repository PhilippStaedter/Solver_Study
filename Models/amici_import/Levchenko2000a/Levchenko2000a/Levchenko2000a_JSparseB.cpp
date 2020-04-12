#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_Levchenko2000a(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 1.0*dwdx0 + 1.0*dwdx1;
    JSparseB[1] = 1.0*dwdx35;
    JSparseB[2] = 1.0*dwdx49;
    JSparseB[3] = -1.0*dwdx148 - 1.0*dwdx149;
    JSparseB[4] = -1.0*dwdx150 - 1.0*dwdx151;
    JSparseB[5] = 1.0*dwdx2 + 1.0*dwdx3;
    JSparseB[6] = 1.0*dwdx75;
    JSparseB[7] = 1.0*dwdx89;
    JSparseB[8] = -1.0*dwdx152 - 1.0*dwdx153;
    JSparseB[9] = -1.0*dwdx154 - 1.0*dwdx155;
    JSparseB[10] = 1.0*dwdx10 + 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15 + 1.0*dwdx16 + 1.0*dwdx17 + 1.0*dwdx18 + 1.0*dwdx19 + 1.0*dwdx20 + 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8 + 1.0*dwdx9;
    JSparseB[11] = 1.0*dwdx104;
    JSparseB[12] = -1.0*dwdx156 - 1.0*dwdx157;
    JSparseB[13] = 1.0*dwdx175;
    JSparseB[14] = 1.0*dwdx194;
    JSparseB[15] = 1.0*dwdx212;
    JSparseB[16] = 1.0*dwdx230;
    JSparseB[17] = 1.0*dwdx247;
    JSparseB[18] = 1.0*dwdx260;
    JSparseB[19] = 1.0*dwdx272;
    JSparseB[20] = 1.0*dwdx286;
    JSparseB[21] = 1.0*dwdx302;
    JSparseB[22] = 1.0*dwdx315;
    JSparseB[23] = 1.0*dwdx327;
    JSparseB[24] = 1.0*dwdx341;
    JSparseB[25] = 1.0*dwdx357;
    JSparseB[26] = 1.0*dwdx370;
    JSparseB[27] = 1.0*dwdx382;
    JSparseB[28] = 1.0*dwdx394;
    JSparseB[29] = -1.0*dwdx398 - 1.0*dwdx399;
    JSparseB[30] = -1.0*dwdx400 - 1.0*dwdx401;
    JSparseB[31] = -1.0*dwdx402 - 1.0*dwdx403;
    JSparseB[32] = -1.0*dwdx404 - 1.0*dwdx405;
    JSparseB[33] = -1.0*dwdx406 - 1.0*dwdx407;
    JSparseB[34] = -1.0*dwdx408 - 1.0*dwdx409;
    JSparseB[35] = -1.0*dwdx410 - 1.0*dwdx411;
    JSparseB[36] = -1.0*dwdx412 - 1.0*dwdx413;
    JSparseB[37] = -1.0*dwdx414 - 1.0*dwdx415;
    JSparseB[38] = -1.0*dwdx416 - 1.0*dwdx417;
    JSparseB[39] = -1.0*dwdx418 - 1.0*dwdx419;
    JSparseB[40] = -1.0*dwdx420 - 1.0*dwdx421;
    JSparseB[41] = -1.0*dwdx422 - 1.0*dwdx423;
    JSparseB[42] = -1.0*dwdx424 - 1.0*dwdx425;
    JSparseB[43] = -1.0*dwdx426 - 1.0*dwdx427;
    JSparseB[44] = -1.0*dwdx428 - 1.0*dwdx429;
    JSparseB[45] = 1.0*dwdx21;
    JSparseB[46] = 1.0*dwdx121;
    JSparseB[47] = -1.0*dwdx158 - 1.0*dwdx159;
    JSparseB[48] = 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24 + 1.0*dwdx25 + 1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28 + 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31 + 1.0*dwdx32 + 1.0*dwdx33 + 1.0*dwdx34;
    JSparseB[49] = 1.0*dwdx90;
    JSparseB[50] = -1.0*dwdx140;
    JSparseB[51] = -1.0*dwdx149;
    JSparseB[52] = 1.0*dwdx160;
    JSparseB[53] = 1.0*dwdx168;
    JSparseB[54] = 1.0*dwdx176;
    JSparseB[55] = 1.0*dwdx183;
    JSparseB[56] = 1.0*dwdx189;
    JSparseB[57] = 1.0*dwdx195;
    JSparseB[58] = 1.0*dwdx201;
    JSparseB[59] = 1.0*dwdx207;
    JSparseB[60] = 1.0*dwdx213;
    JSparseB[61] = 1.0*dwdx219;
    JSparseB[62] = 1.0*dwdx225;
    JSparseB[63] = 1.0*dwdx231;
    JSparseB[64] = -1.0*dwdx236;
    JSparseB[65] = -1.0*dwdx242;
    JSparseB[66] = -1.0*dwdx248;
    JSparseB[67] = -1.0*dwdx253;
    JSparseB[68] = -1.0*dwdx257;
    JSparseB[69] = -1.0*dwdx261;
    JSparseB[70] = -1.0*dwdx265;
    JSparseB[71] = -1.0*dwdx269;
    JSparseB[72] = -1.0*dwdx273;
    JSparseB[73] = -1.0*dwdx277;
    JSparseB[74] = -1.0*dwdx282;
    JSparseB[75] = -1.0*dwdx287;
    JSparseB[76] = 1.0*dwdx0;
    JSparseB[77] = 1.0*dwdx35 + 1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx38 + 1.0*dwdx39 + 1.0*dwdx40 + 1.0*dwdx41 + 1.0*dwdx42 + 1.0*dwdx43 + 1.0*dwdx44 + 1.0*dwdx45 + 1.0*dwdx46 + 1.0*dwdx47 + 1.0*dwdx48;
    JSparseB[78] = 1.0*dwdx91;
    JSparseB[79] = -1.0*dwdx141;
    JSparseB[80] = -1.0*dwdx142;
    JSparseB[81] = -1.0*dwdx148;
    JSparseB[82] = -1.0*dwdx151;
    JSparseB[83] = 1.0*dwdx161;
    JSparseB[84] = 1.0*dwdx169;
    JSparseB[85] = 1.0*dwdx177;
    JSparseB[86] = 1.0*dwdx184;
    JSparseB[87] = 1.0*dwdx190;
    JSparseB[88] = 1.0*dwdx196;
    JSparseB[89] = 1.0*dwdx202;
    JSparseB[90] = 1.0*dwdx208;
    JSparseB[91] = 1.0*dwdx214;
    JSparseB[92] = 1.0*dwdx220;
    JSparseB[93] = 1.0*dwdx226;
    JSparseB[94] = 1.0*dwdx232;
    JSparseB[95] = -1.0*dwdx291;
    JSparseB[96] = -1.0*dwdx297;
    JSparseB[97] = -1.0*dwdx303;
    JSparseB[98] = -1.0*dwdx308;
    JSparseB[99] = -1.0*dwdx312;
    JSparseB[100] = -1.0*dwdx316;
    JSparseB[101] = -1.0*dwdx320;
    JSparseB[102] = -1.0*dwdx324;
    JSparseB[103] = -1.0*dwdx328;
    JSparseB[104] = -1.0*dwdx332;
    JSparseB[105] = -1.0*dwdx337;
    JSparseB[106] = -1.0*dwdx342;
    JSparseB[107] = 1.0*dwdx1;
    JSparseB[108] = 1.0*dwdx49 + 1.0*dwdx50 + 1.0*dwdx51 + 1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54 + 1.0*dwdx55 + 1.0*dwdx56 + 1.0*dwdx57 + 1.0*dwdx58 + 1.0*dwdx59 + 1.0*dwdx60 + 1.0*dwdx61;
    JSparseB[109] = -1.0*dwdx143;
    JSparseB[110] = -1.0*dwdx150;
    JSparseB[111] = 1.0*dwdx162;
    JSparseB[112] = 1.0*dwdx170;
    JSparseB[113] = 1.0*dwdx178;
    JSparseB[114] = 1.0*dwdx185;
    JSparseB[115] = 1.0*dwdx191;
    JSparseB[116] = 1.0*dwdx197;
    JSparseB[117] = 1.0*dwdx203;
    JSparseB[118] = 1.0*dwdx209;
    JSparseB[119] = 1.0*dwdx215;
    JSparseB[120] = 1.0*dwdx221;
    JSparseB[121] = 1.0*dwdx227;
    JSparseB[122] = 1.0*dwdx233;
    JSparseB[123] = -1.0*dwdx346;
    JSparseB[124] = -1.0*dwdx352;
    JSparseB[125] = -1.0*dwdx358;
    JSparseB[126] = -1.0*dwdx363;
    JSparseB[127] = -1.0*dwdx367;
    JSparseB[128] = -1.0*dwdx371;
    JSparseB[129] = -1.0*dwdx375;
    JSparseB[130] = -1.0*dwdx379;
    JSparseB[131] = -1.0*dwdx383;
    JSparseB[132] = -1.0*dwdx387;
    JSparseB[133] = -1.0*dwdx391;
    JSparseB[134] = -1.0*dwdx395;
    JSparseB[135] = 1.0*dwdx62 + 1.0*dwdx63 + 1.0*dwdx64 + 1.0*dwdx65 + 1.0*dwdx66 + 1.0*dwdx67 + 1.0*dwdx68 + 1.0*dwdx69 + 1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72 + 1.0*dwdx73 + 1.0*dwdx74;
    JSparseB[136] = 1.0*dwdx122;
    JSparseB[137] = -1.0*dwdx144;
    JSparseB[138] = -1.0*dwdx153;
    JSparseB[139] = 1.0*dwdx163;
    JSparseB[140] = 1.0*dwdx171;
    JSparseB[141] = 1.0*dwdx179;
    JSparseB[142] = -1.0*dwdx186;
    JSparseB[143] = -1.0*dwdx192;
    JSparseB[144] = -1.0*dwdx198;
    JSparseB[145] = 1.0*dwdx237;
    JSparseB[146] = 1.0*dwdx243;
    JSparseB[147] = 1.0*dwdx249;
    JSparseB[148] = -1.0*dwdx254;
    JSparseB[149] = -1.0*dwdx258;
    JSparseB[150] = -1.0*dwdx262;
    JSparseB[151] = 1.0*dwdx292;
    JSparseB[152] = 1.0*dwdx298;
    JSparseB[153] = 1.0*dwdx304;
    JSparseB[154] = -1.0*dwdx309;
    JSparseB[155] = -1.0*dwdx313;
    JSparseB[156] = -1.0*dwdx317;
    JSparseB[157] = 1.0*dwdx347;
    JSparseB[158] = 1.0*dwdx353;
    JSparseB[159] = 1.0*dwdx359;
    JSparseB[160] = -1.0*dwdx364;
    JSparseB[161] = -1.0*dwdx368;
    JSparseB[162] = -1.0*dwdx372;
    JSparseB[163] = 1.0*dwdx2;
    JSparseB[164] = 1.0*dwdx75 + 1.0*dwdx76 + 1.0*dwdx77 + 1.0*dwdx78 + 1.0*dwdx79 + 1.0*dwdx80 + 1.0*dwdx81 + 1.0*dwdx82 + 1.0*dwdx83 + 1.0*dwdx84 + 1.0*dwdx85 + 1.0*dwdx86 + 1.0*dwdx87 + 1.0*dwdx88;
    JSparseB[165] = 1.0*dwdx123;
    JSparseB[166] = -1.0*dwdx145;
    JSparseB[167] = -1.0*dwdx146;
    JSparseB[168] = -1.0*dwdx152;
    JSparseB[169] = -1.0*dwdx155;
    JSparseB[170] = 1.0*dwdx164;
    JSparseB[171] = 1.0*dwdx172;
    JSparseB[172] = 1.0*dwdx180;
    JSparseB[173] = -1.0*dwdx204;
    JSparseB[174] = -1.0*dwdx210;
    JSparseB[175] = -1.0*dwdx216;
    JSparseB[176] = 1.0*dwdx238;
    JSparseB[177] = 1.0*dwdx244;
    JSparseB[178] = 1.0*dwdx250;
    JSparseB[179] = -1.0*dwdx266;
    JSparseB[180] = -1.0*dwdx270;
    JSparseB[181] = -1.0*dwdx274;
    JSparseB[182] = 1.0*dwdx293;
    JSparseB[183] = 1.0*dwdx299;
    JSparseB[184] = 1.0*dwdx305;
    JSparseB[185] = -1.0*dwdx321;
    JSparseB[186] = -1.0*dwdx325;
    JSparseB[187] = -1.0*dwdx329;
    JSparseB[188] = 1.0*dwdx348;
    JSparseB[189] = 1.0*dwdx354;
    JSparseB[190] = 1.0*dwdx360;
    JSparseB[191] = -1.0*dwdx376;
    JSparseB[192] = -1.0*dwdx380;
    JSparseB[193] = -1.0*dwdx384;
    JSparseB[194] = 1.0*dwdx3;
    JSparseB[195] = 1.0*dwdx22;
    JSparseB[196] = 1.0*dwdx36;
    JSparseB[197] = 1.0*dwdx100 + 1.0*dwdx101 + 1.0*dwdx102 + 1.0*dwdx103 + 1.0*dwdx89 + 1.0*dwdx90 + 1.0*dwdx91 + 1.0*dwdx92 + 1.0*dwdx93 + 1.0*dwdx94 + 1.0*dwdx95 + 1.0*dwdx96 + 1.0*dwdx97 + 1.0*dwdx98 + 1.0*dwdx99;
    JSparseB[198] = -1.0*dwdx140 - 1.0*dwdx141;
    JSparseB[199] = -1.0*dwdx142 - 1.0*dwdx143;
    JSparseB[200] = -1.0*dwdx147;
    JSparseB[201] = -1.0*dwdx154;
    JSparseB[202] = 1.0*dwdx165;
    JSparseB[203] = 1.0*dwdx173;
    JSparseB[204] = 1.0*dwdx181;
    JSparseB[205] = -1.0*dwdx222;
    JSparseB[206] = -1.0*dwdx228;
    JSparseB[207] = -1.0*dwdx234;
    JSparseB[208] = 1.0*dwdx239;
    JSparseB[209] = 1.0*dwdx245;
    JSparseB[210] = 1.0*dwdx251;
    JSparseB[211] = -1.0*dwdx278;
    JSparseB[212] = -1.0*dwdx283;
    JSparseB[213] = -1.0*dwdx288;
    JSparseB[214] = 1.0*dwdx294;
    JSparseB[215] = 1.0*dwdx300;
    JSparseB[216] = 1.0*dwdx306;
    JSparseB[217] = -1.0*dwdx333;
    JSparseB[218] = -1.0*dwdx338;
    JSparseB[219] = -1.0*dwdx343;
    JSparseB[220] = 1.0*dwdx349;
    JSparseB[221] = 1.0*dwdx355;
    JSparseB[222] = 1.0*dwdx361;
    JSparseB[223] = -1.0*dwdx388;
    JSparseB[224] = -1.0*dwdx392;
    JSparseB[225] = -1.0*dwdx396;
    JSparseB[226] = 1.0*dwdx4;
    JSparseB[227] = 1.0*dwdx104 + 1.0*dwdx105 + 1.0*dwdx106 + 1.0*dwdx107 + 1.0*dwdx108 + 1.0*dwdx109 + 1.0*dwdx110 + 1.0*dwdx111 + 1.0*dwdx112 + 1.0*dwdx113 + 1.0*dwdx114 + 1.0*dwdx115 + 1.0*dwdx116 + 1.0*dwdx117 + 1.0*dwdx118 + 1.0*dwdx119 + 1.0*dwdx120;
    JSparseB[228] = -1.0*dwdx156;
    JSparseB[229] = -1.0*dwdx159;
    JSparseB[230] = 1.0*dwdx166;
    JSparseB[231] = -1.0*dwdx174;
    JSparseB[232] = 1.0*dwdx187;
    JSparseB[233] = -1.0*dwdx193;
    JSparseB[234] = 1.0*dwdx205;
    JSparseB[235] = -1.0*dwdx211;
    JSparseB[236] = 1.0*dwdx223;
    JSparseB[237] = -1.0*dwdx229;
    JSparseB[238] = 1.0*dwdx240;
    JSparseB[239] = -1.0*dwdx246;
    JSparseB[240] = 1.0*dwdx255;
    JSparseB[241] = -1.0*dwdx259;
    JSparseB[242] = 1.0*dwdx267;
    JSparseB[243] = -1.0*dwdx271;
    JSparseB[244] = 1.0*dwdx279;
    JSparseB[245] = -1.0*dwdx284;
    JSparseB[246] = 1.0*dwdx295;
    JSparseB[247] = -1.0*dwdx301;
    JSparseB[248] = 1.0*dwdx310;
    JSparseB[249] = -1.0*dwdx314;
    JSparseB[250] = 1.0*dwdx322;
    JSparseB[251] = -1.0*dwdx326;
    JSparseB[252] = 1.0*dwdx334;
    JSparseB[253] = -1.0*dwdx339;
    JSparseB[254] = 1.0*dwdx350;
    JSparseB[255] = -1.0*dwdx356;
    JSparseB[256] = 1.0*dwdx365;
    JSparseB[257] = -1.0*dwdx369;
    JSparseB[258] = 1.0*dwdx377;
    JSparseB[259] = -1.0*dwdx381;
    JSparseB[260] = 1.0*dwdx389;
    JSparseB[261] = -1.0*dwdx393;
    JSparseB[262] = 1.0*dwdx21;
    JSparseB[263] = 1.0*dwdx62;
    JSparseB[264] = 1.0*dwdx76;
    JSparseB[265] = 1.0*dwdx121 + 1.0*dwdx122 + 1.0*dwdx123 + 1.0*dwdx124 + 1.0*dwdx125 + 1.0*dwdx126 + 1.0*dwdx127 + 1.0*dwdx128 + 1.0*dwdx129 + 1.0*dwdx130 + 1.0*dwdx131 + 1.0*dwdx132 + 1.0*dwdx133 + 1.0*dwdx134 + 1.0*dwdx135 + 1.0*dwdx136 + 1.0*dwdx137 + 1.0*dwdx138 + 1.0*dwdx139;
    JSparseB[266] = -1.0*dwdx144 - 1.0*dwdx145;
    JSparseB[267] = -1.0*dwdx146 - 1.0*dwdx147;
    JSparseB[268] = -1.0*dwdx157;
    JSparseB[269] = -1.0*dwdx158;
    JSparseB[270] = 1.0*dwdx167;
    JSparseB[271] = -1.0*dwdx182;
    JSparseB[272] = 1.0*dwdx188;
    JSparseB[273] = -1.0*dwdx199;
    JSparseB[274] = 1.0*dwdx206;
    JSparseB[275] = -1.0*dwdx217;
    JSparseB[276] = 1.0*dwdx224;
    JSparseB[277] = -1.0*dwdx235;
    JSparseB[278] = 1.0*dwdx241;
    JSparseB[279] = -1.0*dwdx252;
    JSparseB[280] = 1.0*dwdx256;
    JSparseB[281] = -1.0*dwdx263;
    JSparseB[282] = 1.0*dwdx268;
    JSparseB[283] = -1.0*dwdx275;
    JSparseB[284] = 1.0*dwdx280;
    JSparseB[285] = -1.0*dwdx289;
    JSparseB[286] = 1.0*dwdx296;
    JSparseB[287] = -1.0*dwdx307;
    JSparseB[288] = 1.0*dwdx311;
    JSparseB[289] = -1.0*dwdx318;
    JSparseB[290] = 1.0*dwdx323;
    JSparseB[291] = -1.0*dwdx330;
    JSparseB[292] = 1.0*dwdx335;
    JSparseB[293] = -1.0*dwdx344;
    JSparseB[294] = 1.0*dwdx351;
    JSparseB[295] = -1.0*dwdx362;
    JSparseB[296] = 1.0*dwdx366;
    JSparseB[297] = -1.0*dwdx373;
    JSparseB[298] = 1.0*dwdx378;
    JSparseB[299] = -1.0*dwdx385;
    JSparseB[300] = 1.0*dwdx390;
    JSparseB[301] = -1.0*dwdx397;
    JSparseB[302] = -1.0*dwdx22;
    JSparseB[303] = -1.0*dwdx90;
    JSparseB[304] = 1.0*dwdx140 + 1.0*dwdx141;
    JSparseB[305] = -1.0*dwdx36;
    JSparseB[306] = -1.0*dwdx91;
    JSparseB[307] = 1.0*dwdx142 + 1.0*dwdx143;
    JSparseB[308] = -1.0*dwdx62;
    JSparseB[309] = -1.0*dwdx122;
    JSparseB[310] = 1.0*dwdx144 + 1.0*dwdx145;
    JSparseB[311] = -1.0*dwdx76;
    JSparseB[312] = -1.0*dwdx123;
    JSparseB[313] = 1.0*dwdx146 + 1.0*dwdx147;
    JSparseB[314] = -1.0*dwdx0;
    JSparseB[315] = -1.0*dwdx35;
    JSparseB[316] = 1.0*dwdx148 + 1.0*dwdx149;
    JSparseB[317] = -1.0*dwdx1;
    JSparseB[318] = -1.0*dwdx49;
    JSparseB[319] = 1.0*dwdx150 + 1.0*dwdx151;
    JSparseB[320] = -1.0*dwdx2;
    JSparseB[321] = -1.0*dwdx75;
    JSparseB[322] = 1.0*dwdx152 + 1.0*dwdx153;
    JSparseB[323] = -1.0*dwdx3;
    JSparseB[324] = -1.0*dwdx89;
    JSparseB[325] = 1.0*dwdx154 + 1.0*dwdx155;
    JSparseB[326] = -1.0*dwdx4;
    JSparseB[327] = -1.0*dwdx104;
    JSparseB[328] = 1.0*dwdx156 + 1.0*dwdx157;
    JSparseB[329] = -1.0*dwdx21;
    JSparseB[330] = -1.0*dwdx121;
    JSparseB[331] = 1.0*dwdx158 + 1.0*dwdx159;
    JSparseB[332] = 1.0*dwdx23;
    JSparseB[333] = 1.0*dwdx37;
    JSparseB[334] = 1.0*dwdx50;
    JSparseB[335] = 1.0*dwdx63;
    JSparseB[336] = 1.0*dwdx77;
    JSparseB[337] = 1.0*dwdx92;
    JSparseB[338] = 1.0*dwdx105;
    JSparseB[339] = 1.0*dwdx124;
    JSparseB[340] = 1.0*dwdx160 + 1.0*dwdx161 + 1.0*dwdx162 + 1.0*dwdx163 + 1.0*dwdx164 + 1.0*dwdx165 + 1.0*dwdx166 + 1.0*dwdx167;
    JSparseB[341] = -1.0*dwdx174;
    JSparseB[342] = -1.0*dwdx182;
    JSparseB[343] = -1.0*dwdx186;
    JSparseB[344] = -1.0*dwdx204;
    JSparseB[345] = -1.0*dwdx222;
    JSparseB[346] = -1.0*dwdx236;
    JSparseB[347] = -1.0*dwdx291;
    JSparseB[348] = -1.0*dwdx346;
    JSparseB[349] = 1.0*dwdx5;
    JSparseB[350] = 1.0*dwdx24;
    JSparseB[351] = 1.0*dwdx38;
    JSparseB[352] = 1.0*dwdx51;
    JSparseB[353] = 1.0*dwdx64;
    JSparseB[354] = 1.0*dwdx78;
    JSparseB[355] = 1.0*dwdx93;
    JSparseB[356] = -1.0*dwdx105;
    JSparseB[357] = -1.0*dwdx166;
    JSparseB[358] = 1.0*dwdx168 + 1.0*dwdx169 + 1.0*dwdx170 + 1.0*dwdx171 + 1.0*dwdx172 + 1.0*dwdx173 + 1.0*dwdx174 + 1.0*dwdx175;
    JSparseB[359] = -1.0*dwdx192;
    JSparseB[360] = -1.0*dwdx210;
    JSparseB[361] = -1.0*dwdx228;
    JSparseB[362] = -1.0*dwdx242;
    JSparseB[363] = -1.0*dwdx297;
    JSparseB[364] = -1.0*dwdx352;
    JSparseB[365] = -1.0*dwdx398;
    JSparseB[366] = 1.0*dwdx25;
    JSparseB[367] = 1.0*dwdx39;
    JSparseB[368] = 1.0*dwdx52;
    JSparseB[369] = 1.0*dwdx65;
    JSparseB[370] = 1.0*dwdx79;
    JSparseB[371] = 1.0*dwdx94;
    JSparseB[372] = -1.0*dwdx124;
    JSparseB[373] = -1.0*dwdx167;
    JSparseB[374] = 1.0*dwdx176 + 1.0*dwdx177 + 1.0*dwdx178 + 1.0*dwdx179 + 1.0*dwdx180 + 1.0*dwdx181 + 1.0*dwdx182;
    JSparseB[375] = -1.0*dwdx198;
    JSparseB[376] = -1.0*dwdx216;
    JSparseB[377] = -1.0*dwdx234;
    JSparseB[378] = -1.0*dwdx248;
    JSparseB[379] = -1.0*dwdx303;
    JSparseB[380] = -1.0*dwdx358;
    JSparseB[381] = -1.0*dwdx399;
    JSparseB[382] = 1.0*dwdx26;
    JSparseB[383] = 1.0*dwdx40;
    JSparseB[384] = 1.0*dwdx53;
    JSparseB[385] = -1.0*dwdx63;
    JSparseB[386] = 1.0*dwdx106;
    JSparseB[387] = 1.0*dwdx125;
    JSparseB[388] = -1.0*dwdx163;
    JSparseB[389] = 1.0*dwdx183 + 1.0*dwdx184 + 1.0*dwdx185 + 1.0*dwdx186 + 1.0*dwdx187 + 1.0*dwdx188;
    JSparseB[390] = -1.0*dwdx193;
    JSparseB[391] = -1.0*dwdx199;
    JSparseB[392] = -1.0*dwdx253;
    JSparseB[393] = -1.0*dwdx308;
    JSparseB[394] = -1.0*dwdx363;
    JSparseB[395] = 1.0*dwdx6;
    JSparseB[396] = 1.0*dwdx27;
    JSparseB[397] = 1.0*dwdx41;
    JSparseB[398] = 1.0*dwdx54;
    JSparseB[399] = -1.0*dwdx64;
    JSparseB[400] = -1.0*dwdx106;
    JSparseB[401] = -1.0*dwdx171;
    JSparseB[402] = -1.0*dwdx187;
    JSparseB[403] = 1.0*dwdx189 + 1.0*dwdx190 + 1.0*dwdx191 + 1.0*dwdx192 + 1.0*dwdx193 + 1.0*dwdx194;
    JSparseB[404] = -1.0*dwdx257;
    JSparseB[405] = -1.0*dwdx312;
    JSparseB[406] = -1.0*dwdx367;
    JSparseB[407] = -1.0*dwdx400;
    JSparseB[408] = 1.0*dwdx28;
    JSparseB[409] = 1.0*dwdx42;
    JSparseB[410] = 1.0*dwdx55;
    JSparseB[411] = -1.0*dwdx65;
    JSparseB[412] = -1.0*dwdx125;
    JSparseB[413] = -1.0*dwdx179;
    JSparseB[414] = -1.0*dwdx188;
    JSparseB[415] = 1.0*dwdx195 + 1.0*dwdx196 + 1.0*dwdx197 + 1.0*dwdx198 + 1.0*dwdx199 + 1.0*dwdx200;
    JSparseB[416] = -1.0*dwdx261;
    JSparseB[417] = -1.0*dwdx316;
    JSparseB[418] = -1.0*dwdx371;
    JSparseB[419] = -1.0*dwdx401;
    JSparseB[420] = 1.0*dwdx29;
    JSparseB[421] = 1.0*dwdx43;
    JSparseB[422] = 1.0*dwdx56;
    JSparseB[423] = -1.0*dwdx77;
    JSparseB[424] = 1.0*dwdx107;
    JSparseB[425] = 1.0*dwdx126;
    JSparseB[426] = -1.0*dwdx164;
    JSparseB[427] = 1.0*dwdx201 + 1.0*dwdx202 + 1.0*dwdx203 + 1.0*dwdx204 + 1.0*dwdx205 + 1.0*dwdx206;
    JSparseB[428] = -1.0*dwdx211;
    JSparseB[429] = -1.0*dwdx217;
    JSparseB[430] = -1.0*dwdx265;
    JSparseB[431] = -1.0*dwdx320;
    JSparseB[432] = -1.0*dwdx375;
    JSparseB[433] = 1.0*dwdx7;
    JSparseB[434] = 1.0*dwdx30;
    JSparseB[435] = 1.0*dwdx44;
    JSparseB[436] = 1.0*dwdx57;
    JSparseB[437] = -1.0*dwdx78;
    JSparseB[438] = -1.0*dwdx107;
    JSparseB[439] = -1.0*dwdx172;
    JSparseB[440] = -1.0*dwdx205;
    JSparseB[441] = 1.0*dwdx207 + 1.0*dwdx208 + 1.0*dwdx209 + 1.0*dwdx210 + 1.0*dwdx211 + 1.0*dwdx212;
    JSparseB[442] = -1.0*dwdx269;
    JSparseB[443] = -1.0*dwdx324;
    JSparseB[444] = -1.0*dwdx379;
    JSparseB[445] = -1.0*dwdx402;
    JSparseB[446] = 1.0*dwdx31;
    JSparseB[447] = 1.0*dwdx45;
    JSparseB[448] = 1.0*dwdx58;
    JSparseB[449] = -1.0*dwdx79;
    JSparseB[450] = -1.0*dwdx126;
    JSparseB[451] = -1.0*dwdx180;
    JSparseB[452] = -1.0*dwdx200;
    JSparseB[453] = -1.0*dwdx206;
    JSparseB[454] = 1.0*dwdx213 + 1.0*dwdx214 + 1.0*dwdx215 + 1.0*dwdx216 + 1.0*dwdx217 + 1.0*dwdx218;
    JSparseB[455] = -1.0*dwdx273;
    JSparseB[456] = -1.0*dwdx328;
    JSparseB[457] = -1.0*dwdx383;
    JSparseB[458] = -1.0*dwdx403;
    JSparseB[459] = 1.0*dwdx32;
    JSparseB[460] = 1.0*dwdx46;
    JSparseB[461] = 1.0*dwdx59;
    JSparseB[462] = -1.0*dwdx92;
    JSparseB[463] = 1.0*dwdx108;
    JSparseB[464] = 1.0*dwdx127;
    JSparseB[465] = -1.0*dwdx165;
    JSparseB[466] = 1.0*dwdx219 + 1.0*dwdx220 + 1.0*dwdx221 + 1.0*dwdx222 + 1.0*dwdx223 + 1.0*dwdx224;
    JSparseB[467] = -1.0*dwdx229;
    JSparseB[468] = -1.0*dwdx235;
    JSparseB[469] = -1.0*dwdx277;
    JSparseB[470] = -1.0*dwdx332;
    JSparseB[471] = -1.0*dwdx387;
    JSparseB[472] = 1.0*dwdx8;
    JSparseB[473] = 1.0*dwdx33;
    JSparseB[474] = 1.0*dwdx47;
    JSparseB[475] = 1.0*dwdx60;
    JSparseB[476] = -1.0*dwdx93;
    JSparseB[477] = -1.0*dwdx108;
    JSparseB[478] = -1.0*dwdx173;
    JSparseB[479] = -1.0*dwdx223;
    JSparseB[480] = 1.0*dwdx225 + 1.0*dwdx226 + 1.0*dwdx227 + 1.0*dwdx228 + 1.0*dwdx229 + 1.0*dwdx230;
    JSparseB[481] = -1.0*dwdx282;
    JSparseB[482] = -1.0*dwdx337;
    JSparseB[483] = -1.0*dwdx391;
    JSparseB[484] = -1.0*dwdx404;
    JSparseB[485] = 1.0*dwdx34;
    JSparseB[486] = 1.0*dwdx48;
    JSparseB[487] = 1.0*dwdx61;
    JSparseB[488] = -1.0*dwdx94;
    JSparseB[489] = -1.0*dwdx127;
    JSparseB[490] = -1.0*dwdx181;
    JSparseB[491] = -1.0*dwdx218;
    JSparseB[492] = -1.0*dwdx224;
    JSparseB[493] = 1.0*dwdx231 + 1.0*dwdx232 + 1.0*dwdx233 + 1.0*dwdx234 + 1.0*dwdx235;
    JSparseB[494] = -1.0*dwdx287;
    JSparseB[495] = -1.0*dwdx342;
    JSparseB[496] = -1.0*dwdx395;
    JSparseB[497] = -1.0*dwdx405;
    JSparseB[498] = -1.0*dwdx23;
    JSparseB[499] = 1.0*dwdx66;
    JSparseB[500] = 1.0*dwdx80;
    JSparseB[501] = 1.0*dwdx95;
    JSparseB[502] = 1.0*dwdx109;
    JSparseB[503] = 1.0*dwdx128;
    JSparseB[504] = -1.0*dwdx160;
    JSparseB[505] = 1.0*dwdx236 + 1.0*dwdx237 + 1.0*dwdx238 + 1.0*dwdx239 + 1.0*dwdx240 + 1.0*dwdx241;
    JSparseB[506] = -1.0*dwdx246;
    JSparseB[507] = -1.0*dwdx252;
    JSparseB[508] = -1.0*dwdx254;
    JSparseB[509] = -1.0*dwdx266;
    JSparseB[510] = -1.0*dwdx278;
    JSparseB[511] = 1.0*dwdx9;
    JSparseB[512] = -1.0*dwdx24;
    JSparseB[513] = 1.0*dwdx67;
    JSparseB[514] = 1.0*dwdx81;
    JSparseB[515] = 1.0*dwdx96;
    JSparseB[516] = -1.0*dwdx109;
    JSparseB[517] = -1.0*dwdx168;
    JSparseB[518] = -1.0*dwdx240;
    JSparseB[519] = 1.0*dwdx242 + 1.0*dwdx243 + 1.0*dwdx244 + 1.0*dwdx245 + 1.0*dwdx246 + 1.0*dwdx247;
    JSparseB[520] = -1.0*dwdx258;
    JSparseB[521] = -1.0*dwdx270;
    JSparseB[522] = -1.0*dwdx283;
    JSparseB[523] = -1.0*dwdx406;
    JSparseB[524] = -1.0*dwdx25;
    JSparseB[525] = 1.0*dwdx68;
    JSparseB[526] = 1.0*dwdx82;
    JSparseB[527] = 1.0*dwdx97;
    JSparseB[528] = -1.0*dwdx128;
    JSparseB[529] = -1.0*dwdx176;
    JSparseB[530] = -1.0*dwdx241;
    JSparseB[531] = 1.0*dwdx248 + 1.0*dwdx249 + 1.0*dwdx250 + 1.0*dwdx251 + 1.0*dwdx252;
    JSparseB[532] = -1.0*dwdx262;
    JSparseB[533] = -1.0*dwdx274;
    JSparseB[534] = -1.0*dwdx288;
    JSparseB[535] = -1.0*dwdx407;
    JSparseB[536] = -1.0*dwdx26;
    JSparseB[537] = -1.0*dwdx66;
    JSparseB[538] = 1.0*dwdx110;
    JSparseB[539] = 1.0*dwdx129;
    JSparseB[540] = -1.0*dwdx183;
    JSparseB[541] = -1.0*dwdx237;
    JSparseB[542] = 1.0*dwdx253 + 1.0*dwdx254 + 1.0*dwdx255 + 1.0*dwdx256;
    JSparseB[543] = -1.0*dwdx259;
    JSparseB[544] = -1.0*dwdx263;
    JSparseB[545] = 1.0*dwdx10;
    JSparseB[546] = -1.0*dwdx27;
    JSparseB[547] = -1.0*dwdx67;
    JSparseB[548] = -1.0*dwdx110;
    JSparseB[549] = -1.0*dwdx189;
    JSparseB[550] = -1.0*dwdx243;
    JSparseB[551] = -1.0*dwdx255;
    JSparseB[552] = 1.0*dwdx257 + 1.0*dwdx258 + 1.0*dwdx259 + 1.0*dwdx260;
    JSparseB[553] = -1.0*dwdx408;
    JSparseB[554] = -1.0*dwdx28;
    JSparseB[555] = -1.0*dwdx68;
    JSparseB[556] = -1.0*dwdx129;
    JSparseB[557] = -1.0*dwdx195;
    JSparseB[558] = -1.0*dwdx249;
    JSparseB[559] = -1.0*dwdx256;
    JSparseB[560] = 1.0*dwdx261 + 1.0*dwdx262 + 1.0*dwdx263 + 1.0*dwdx264;
    JSparseB[561] = -1.0*dwdx409;
    JSparseB[562] = -1.0*dwdx29;
    JSparseB[563] = -1.0*dwdx80;
    JSparseB[564] = 1.0*dwdx111;
    JSparseB[565] = 1.0*dwdx130;
    JSparseB[566] = -1.0*dwdx201;
    JSparseB[567] = -1.0*dwdx238;
    JSparseB[568] = 1.0*dwdx265 + 1.0*dwdx266 + 1.0*dwdx267 + 1.0*dwdx268;
    JSparseB[569] = -1.0*dwdx271;
    JSparseB[570] = -1.0*dwdx275;
    JSparseB[571] = 1.0*dwdx11;
    JSparseB[572] = -1.0*dwdx30;
    JSparseB[573] = -1.0*dwdx81;
    JSparseB[574] = -1.0*dwdx111;
    JSparseB[575] = -1.0*dwdx207;
    JSparseB[576] = -1.0*dwdx244;
    JSparseB[577] = -1.0*dwdx267;
    JSparseB[578] = 1.0*dwdx269 + 1.0*dwdx270 + 1.0*dwdx271 + 1.0*dwdx272;
    JSparseB[579] = -1.0*dwdx410;
    JSparseB[580] = -1.0*dwdx31;
    JSparseB[581] = -1.0*dwdx82;
    JSparseB[582] = -1.0*dwdx130;
    JSparseB[583] = -1.0*dwdx213;
    JSparseB[584] = -1.0*dwdx250;
    JSparseB[585] = -1.0*dwdx264;
    JSparseB[586] = -1.0*dwdx268;
    JSparseB[587] = 1.0*dwdx273 + 1.0*dwdx274 + 1.0*dwdx275 + 1.0*dwdx276;
    JSparseB[588] = -1.0*dwdx411;
    JSparseB[589] = -1.0*dwdx32;
    JSparseB[590] = -1.0*dwdx95;
    JSparseB[591] = 1.0*dwdx112;
    JSparseB[592] = 1.0*dwdx131;
    JSparseB[593] = -1.0*dwdx219;
    JSparseB[594] = -1.0*dwdx239;
    JSparseB[595] = 1.0*dwdx277 + 1.0*dwdx278 + 1.0*dwdx279 + 1.0*dwdx280 + 1.0*dwdx281;
    JSparseB[596] = -1.0*dwdx284;
    JSparseB[597] = -1.0*dwdx289;
    JSparseB[598] = 1.0*dwdx12;
    JSparseB[599] = -1.0*dwdx33;
    JSparseB[600] = -1.0*dwdx96;
    JSparseB[601] = -1.0*dwdx112;
    JSparseB[602] = -1.0*dwdx225;
    JSparseB[603] = -1.0*dwdx245;
    JSparseB[604] = -1.0*dwdx279;
    JSparseB[605] = 1.0*dwdx282 + 1.0*dwdx283 + 1.0*dwdx284 + 1.0*dwdx285 + 1.0*dwdx286;
    JSparseB[606] = -1.0*dwdx412;
    JSparseB[607] = -1.0*dwdx34;
    JSparseB[608] = -1.0*dwdx97;
    JSparseB[609] = -1.0*dwdx131;
    JSparseB[610] = -1.0*dwdx231;
    JSparseB[611] = -1.0*dwdx251;
    JSparseB[612] = -1.0*dwdx276;
    JSparseB[613] = -1.0*dwdx280;
    JSparseB[614] = 1.0*dwdx287 + 1.0*dwdx288 + 1.0*dwdx289 + 1.0*dwdx290;
    JSparseB[615] = -1.0*dwdx413;
    JSparseB[616] = -1.0*dwdx37;
    JSparseB[617] = 1.0*dwdx69;
    JSparseB[618] = 1.0*dwdx83;
    JSparseB[619] = 1.0*dwdx98;
    JSparseB[620] = 1.0*dwdx113;
    JSparseB[621] = 1.0*dwdx132;
    JSparseB[622] = -1.0*dwdx161;
    JSparseB[623] = 1.0*dwdx291 + 1.0*dwdx292 + 1.0*dwdx293 + 1.0*dwdx294 + 1.0*dwdx295 + 1.0*dwdx296;
    JSparseB[624] = -1.0*dwdx301;
    JSparseB[625] = -1.0*dwdx307;
    JSparseB[626] = -1.0*dwdx309;
    JSparseB[627] = -1.0*dwdx321;
    JSparseB[628] = -1.0*dwdx333;
    JSparseB[629] = 1.0*dwdx13;
    JSparseB[630] = -1.0*dwdx38;
    JSparseB[631] = 1.0*dwdx70;
    JSparseB[632] = 1.0*dwdx84;
    JSparseB[633] = 1.0*dwdx99;
    JSparseB[634] = -1.0*dwdx113;
    JSparseB[635] = -1.0*dwdx169;
    JSparseB[636] = -1.0*dwdx295;
    JSparseB[637] = 1.0*dwdx297 + 1.0*dwdx298 + 1.0*dwdx299 + 1.0*dwdx300 + 1.0*dwdx301 + 1.0*dwdx302;
    JSparseB[638] = -1.0*dwdx313;
    JSparseB[639] = -1.0*dwdx325;
    JSparseB[640] = -1.0*dwdx338;
    JSparseB[641] = -1.0*dwdx414;
    JSparseB[642] = -1.0*dwdx39;
    JSparseB[643] = 1.0*dwdx71;
    JSparseB[644] = 1.0*dwdx85;
    JSparseB[645] = 1.0*dwdx100;
    JSparseB[646] = -1.0*dwdx132;
    JSparseB[647] = -1.0*dwdx177;
    JSparseB[648] = -1.0*dwdx296;
    JSparseB[649] = 1.0*dwdx303 + 1.0*dwdx304 + 1.0*dwdx305 + 1.0*dwdx306 + 1.0*dwdx307;
    JSparseB[650] = -1.0*dwdx317;
    JSparseB[651] = -1.0*dwdx329;
    JSparseB[652] = -1.0*dwdx343;
    JSparseB[653] = -1.0*dwdx415;
    JSparseB[654] = -1.0*dwdx40;
    JSparseB[655] = -1.0*dwdx69;
    JSparseB[656] = 1.0*dwdx114;
    JSparseB[657] = 1.0*dwdx133;
    JSparseB[658] = -1.0*dwdx184;
    JSparseB[659] = -1.0*dwdx292;
    JSparseB[660] = 1.0*dwdx308 + 1.0*dwdx309 + 1.0*dwdx310 + 1.0*dwdx311;
    JSparseB[661] = -1.0*dwdx314;
    JSparseB[662] = -1.0*dwdx318;
    JSparseB[663] = 1.0*dwdx14;
    JSparseB[664] = -1.0*dwdx41;
    JSparseB[665] = -1.0*dwdx70;
    JSparseB[666] = -1.0*dwdx114;
    JSparseB[667] = -1.0*dwdx190;
    JSparseB[668] = -1.0*dwdx298;
    JSparseB[669] = -1.0*dwdx310;
    JSparseB[670] = 1.0*dwdx312 + 1.0*dwdx313 + 1.0*dwdx314 + 1.0*dwdx315;
    JSparseB[671] = -1.0*dwdx416;
    JSparseB[672] = -1.0*dwdx42;
    JSparseB[673] = -1.0*dwdx71;
    JSparseB[674] = -1.0*dwdx133;
    JSparseB[675] = -1.0*dwdx196;
    JSparseB[676] = -1.0*dwdx304;
    JSparseB[677] = -1.0*dwdx311;
    JSparseB[678] = 1.0*dwdx316 + 1.0*dwdx317 + 1.0*dwdx318 + 1.0*dwdx319;
    JSparseB[679] = -1.0*dwdx417;
    JSparseB[680] = -1.0*dwdx43;
    JSparseB[681] = -1.0*dwdx83;
    JSparseB[682] = 1.0*dwdx115;
    JSparseB[683] = 1.0*dwdx134;
    JSparseB[684] = -1.0*dwdx202;
    JSparseB[685] = -1.0*dwdx293;
    JSparseB[686] = 1.0*dwdx320 + 1.0*dwdx321 + 1.0*dwdx322 + 1.0*dwdx323;
    JSparseB[687] = -1.0*dwdx326;
    JSparseB[688] = -1.0*dwdx330;
    JSparseB[689] = 1.0*dwdx15;
    JSparseB[690] = -1.0*dwdx44;
    JSparseB[691] = -1.0*dwdx84;
    JSparseB[692] = -1.0*dwdx115;
    JSparseB[693] = -1.0*dwdx208;
    JSparseB[694] = -1.0*dwdx299;
    JSparseB[695] = -1.0*dwdx322;
    JSparseB[696] = 1.0*dwdx324 + 1.0*dwdx325 + 1.0*dwdx326 + 1.0*dwdx327;
    JSparseB[697] = -1.0*dwdx418;
    JSparseB[698] = -1.0*dwdx45;
    JSparseB[699] = -1.0*dwdx85;
    JSparseB[700] = -1.0*dwdx134;
    JSparseB[701] = -1.0*dwdx214;
    JSparseB[702] = -1.0*dwdx305;
    JSparseB[703] = -1.0*dwdx319;
    JSparseB[704] = -1.0*dwdx323;
    JSparseB[705] = 1.0*dwdx328 + 1.0*dwdx329 + 1.0*dwdx330 + 1.0*dwdx331;
    JSparseB[706] = -1.0*dwdx419;
    JSparseB[707] = -1.0*dwdx46;
    JSparseB[708] = -1.0*dwdx98;
    JSparseB[709] = 1.0*dwdx116;
    JSparseB[710] = 1.0*dwdx135;
    JSparseB[711] = -1.0*dwdx220;
    JSparseB[712] = -1.0*dwdx281;
    JSparseB[713] = -1.0*dwdx294;
    JSparseB[714] = 1.0*dwdx332 + 1.0*dwdx333 + 1.0*dwdx334 + 1.0*dwdx335 + 1.0*dwdx336;
    JSparseB[715] = -1.0*dwdx339;
    JSparseB[716] = -1.0*dwdx344;
    JSparseB[717] = 1.0*dwdx16;
    JSparseB[718] = -1.0*dwdx47;
    JSparseB[719] = -1.0*dwdx99;
    JSparseB[720] = -1.0*dwdx116;
    JSparseB[721] = -1.0*dwdx226;
    JSparseB[722] = -1.0*dwdx285;
    JSparseB[723] = -1.0*dwdx300;
    JSparseB[724] = -1.0*dwdx334;
    JSparseB[725] = 1.0*dwdx337 + 1.0*dwdx338 + 1.0*dwdx339 + 1.0*dwdx340 + 1.0*dwdx341;
    JSparseB[726] = -1.0*dwdx420;
    JSparseB[727] = -1.0*dwdx48;
    JSparseB[728] = -1.0*dwdx100;
    JSparseB[729] = -1.0*dwdx135;
    JSparseB[730] = -1.0*dwdx232;
    JSparseB[731] = -1.0*dwdx290;
    JSparseB[732] = -1.0*dwdx306;
    JSparseB[733] = -1.0*dwdx331;
    JSparseB[734] = -1.0*dwdx335;
    JSparseB[735] = 1.0*dwdx342 + 1.0*dwdx343 + 1.0*dwdx344 + 1.0*dwdx345;
    JSparseB[736] = -1.0*dwdx421;
    JSparseB[737] = -1.0*dwdx50;
    JSparseB[738] = 1.0*dwdx72;
    JSparseB[739] = 1.0*dwdx86;
    JSparseB[740] = 1.0*dwdx101;
    JSparseB[741] = 1.0*dwdx117;
    JSparseB[742] = 1.0*dwdx136;
    JSparseB[743] = -1.0*dwdx162;
    JSparseB[744] = 1.0*dwdx346 + 1.0*dwdx347 + 1.0*dwdx348 + 1.0*dwdx349 + 1.0*dwdx350 + 1.0*dwdx351;
    JSparseB[745] = -1.0*dwdx356;
    JSparseB[746] = -1.0*dwdx362;
    JSparseB[747] = -1.0*dwdx364;
    JSparseB[748] = -1.0*dwdx376;
    JSparseB[749] = -1.0*dwdx388;
    JSparseB[750] = 1.0*dwdx17;
    JSparseB[751] = -1.0*dwdx51;
    JSparseB[752] = 1.0*dwdx73;
    JSparseB[753] = 1.0*dwdx87;
    JSparseB[754] = 1.0*dwdx102;
    JSparseB[755] = -1.0*dwdx117;
    JSparseB[756] = -1.0*dwdx170;
    JSparseB[757] = -1.0*dwdx350;
    JSparseB[758] = 1.0*dwdx352 + 1.0*dwdx353 + 1.0*dwdx354 + 1.0*dwdx355 + 1.0*dwdx356 + 1.0*dwdx357;
    JSparseB[759] = -1.0*dwdx368;
    JSparseB[760] = -1.0*dwdx380;
    JSparseB[761] = -1.0*dwdx392;
    JSparseB[762] = -1.0*dwdx422;
    JSparseB[763] = -1.0*dwdx52;
    JSparseB[764] = 1.0*dwdx74;
    JSparseB[765] = 1.0*dwdx88;
    JSparseB[766] = 1.0*dwdx103;
    JSparseB[767] = -1.0*dwdx136;
    JSparseB[768] = -1.0*dwdx178;
    JSparseB[769] = -1.0*dwdx351;
    JSparseB[770] = 1.0*dwdx358 + 1.0*dwdx359 + 1.0*dwdx360 + 1.0*dwdx361 + 1.0*dwdx362;
    JSparseB[771] = -1.0*dwdx372;
    JSparseB[772] = -1.0*dwdx384;
    JSparseB[773] = -1.0*dwdx396;
    JSparseB[774] = -1.0*dwdx423;
    JSparseB[775] = -1.0*dwdx53;
    JSparseB[776] = -1.0*dwdx72;
    JSparseB[777] = 1.0*dwdx118;
    JSparseB[778] = 1.0*dwdx137;
    JSparseB[779] = -1.0*dwdx185;
    JSparseB[780] = -1.0*dwdx347;
    JSparseB[781] = 1.0*dwdx363 + 1.0*dwdx364 + 1.0*dwdx365 + 1.0*dwdx366;
    JSparseB[782] = -1.0*dwdx369;
    JSparseB[783] = -1.0*dwdx373;
    JSparseB[784] = 1.0*dwdx18;
    JSparseB[785] = -1.0*dwdx54;
    JSparseB[786] = -1.0*dwdx73;
    JSparseB[787] = -1.0*dwdx118;
    JSparseB[788] = -1.0*dwdx191;
    JSparseB[789] = -1.0*dwdx353;
    JSparseB[790] = -1.0*dwdx365;
    JSparseB[791] = 1.0*dwdx367 + 1.0*dwdx368 + 1.0*dwdx369 + 1.0*dwdx370;
    JSparseB[792] = -1.0*dwdx424;
    JSparseB[793] = -1.0*dwdx55;
    JSparseB[794] = -1.0*dwdx74;
    JSparseB[795] = -1.0*dwdx137;
    JSparseB[796] = -1.0*dwdx197;
    JSparseB[797] = -1.0*dwdx359;
    JSparseB[798] = -1.0*dwdx366;
    JSparseB[799] = 1.0*dwdx371 + 1.0*dwdx372 + 1.0*dwdx373 + 1.0*dwdx374;
    JSparseB[800] = -1.0*dwdx425;
    JSparseB[801] = -1.0*dwdx56;
    JSparseB[802] = -1.0*dwdx86;
    JSparseB[803] = 1.0*dwdx119;
    JSparseB[804] = 1.0*dwdx138;
    JSparseB[805] = -1.0*dwdx203;
    JSparseB[806] = -1.0*dwdx348;
    JSparseB[807] = 1.0*dwdx375 + 1.0*dwdx376 + 1.0*dwdx377 + 1.0*dwdx378;
    JSparseB[808] = -1.0*dwdx381;
    JSparseB[809] = -1.0*dwdx385;
    JSparseB[810] = 1.0*dwdx19;
    JSparseB[811] = -1.0*dwdx57;
    JSparseB[812] = -1.0*dwdx87;
    JSparseB[813] = -1.0*dwdx119;
    JSparseB[814] = -1.0*dwdx209;
    JSparseB[815] = -1.0*dwdx354;
    JSparseB[816] = -1.0*dwdx377;
    JSparseB[817] = 1.0*dwdx379 + 1.0*dwdx380 + 1.0*dwdx381 + 1.0*dwdx382;
    JSparseB[818] = -1.0*dwdx426;
    JSparseB[819] = -1.0*dwdx58;
    JSparseB[820] = -1.0*dwdx88;
    JSparseB[821] = -1.0*dwdx138;
    JSparseB[822] = -1.0*dwdx215;
    JSparseB[823] = -1.0*dwdx360;
    JSparseB[824] = -1.0*dwdx374;
    JSparseB[825] = -1.0*dwdx378;
    JSparseB[826] = 1.0*dwdx383 + 1.0*dwdx384 + 1.0*dwdx385 + 1.0*dwdx386;
    JSparseB[827] = -1.0*dwdx427;
    JSparseB[828] = -1.0*dwdx59;
    JSparseB[829] = -1.0*dwdx101;
    JSparseB[830] = 1.0*dwdx120;
    JSparseB[831] = 1.0*dwdx139;
    JSparseB[832] = -1.0*dwdx221;
    JSparseB[833] = -1.0*dwdx336;
    JSparseB[834] = -1.0*dwdx349;
    JSparseB[835] = 1.0*dwdx387 + 1.0*dwdx388 + 1.0*dwdx389 + 1.0*dwdx390;
    JSparseB[836] = -1.0*dwdx393;
    JSparseB[837] = -1.0*dwdx397;
    JSparseB[838] = 1.0*dwdx20;
    JSparseB[839] = -1.0*dwdx60;
    JSparseB[840] = -1.0*dwdx102;
    JSparseB[841] = -1.0*dwdx120;
    JSparseB[842] = -1.0*dwdx227;
    JSparseB[843] = -1.0*dwdx340;
    JSparseB[844] = -1.0*dwdx355;
    JSparseB[845] = -1.0*dwdx389;
    JSparseB[846] = 1.0*dwdx391 + 1.0*dwdx392 + 1.0*dwdx393 + 1.0*dwdx394;
    JSparseB[847] = -1.0*dwdx428;
    JSparseB[848] = -1.0*dwdx61;
    JSparseB[849] = -1.0*dwdx103;
    JSparseB[850] = -1.0*dwdx139;
    JSparseB[851] = -1.0*dwdx233;
    JSparseB[852] = -1.0*dwdx345;
    JSparseB[853] = -1.0*dwdx361;
    JSparseB[854] = -1.0*dwdx386;
    JSparseB[855] = -1.0*dwdx390;
    JSparseB[856] = 1.0*dwdx395 + 1.0*dwdx396 + 1.0*dwdx397;
    JSparseB[857] = -1.0*dwdx429;
    JSparseB[858] = -1.0*dwdx5;
    JSparseB[859] = -1.0*dwdx175;
    JSparseB[860] = 1.0*dwdx398 + 1.0*dwdx399;
    JSparseB[861] = -1.0*dwdx6;
    JSparseB[862] = -1.0*dwdx194;
    JSparseB[863] = 1.0*dwdx400 + 1.0*dwdx401;
    JSparseB[864] = -1.0*dwdx7;
    JSparseB[865] = -1.0*dwdx212;
    JSparseB[866] = 1.0*dwdx402 + 1.0*dwdx403;
    JSparseB[867] = -1.0*dwdx8;
    JSparseB[868] = -1.0*dwdx230;
    JSparseB[869] = 1.0*dwdx404 + 1.0*dwdx405;
    JSparseB[870] = -1.0*dwdx9;
    JSparseB[871] = -1.0*dwdx247;
    JSparseB[872] = 1.0*dwdx406 + 1.0*dwdx407;
    JSparseB[873] = -1.0*dwdx10;
    JSparseB[874] = -1.0*dwdx260;
    JSparseB[875] = 1.0*dwdx408 + 1.0*dwdx409;
    JSparseB[876] = -1.0*dwdx11;
    JSparseB[877] = -1.0*dwdx272;
    JSparseB[878] = 1.0*dwdx410 + 1.0*dwdx411;
    JSparseB[879] = -1.0*dwdx12;
    JSparseB[880] = -1.0*dwdx286;
    JSparseB[881] = 1.0*dwdx412 + 1.0*dwdx413;
    JSparseB[882] = -1.0*dwdx13;
    JSparseB[883] = -1.0*dwdx302;
    JSparseB[884] = 1.0*dwdx414 + 1.0*dwdx415;
    JSparseB[885] = -1.0*dwdx14;
    JSparseB[886] = -1.0*dwdx315;
    JSparseB[887] = 1.0*dwdx416 + 1.0*dwdx417;
    JSparseB[888] = -1.0*dwdx15;
    JSparseB[889] = -1.0*dwdx327;
    JSparseB[890] = 1.0*dwdx418 + 1.0*dwdx419;
    JSparseB[891] = -1.0*dwdx16;
    JSparseB[892] = -1.0*dwdx341;
    JSparseB[893] = 1.0*dwdx420 + 1.0*dwdx421;
    JSparseB[894] = -1.0*dwdx17;
    JSparseB[895] = -1.0*dwdx357;
    JSparseB[896] = 1.0*dwdx422 + 1.0*dwdx423;
    JSparseB[897] = -1.0*dwdx18;
    JSparseB[898] = -1.0*dwdx370;
    JSparseB[899] = 1.0*dwdx424 + 1.0*dwdx425;
    JSparseB[900] = -1.0*dwdx19;
    JSparseB[901] = -1.0*dwdx382;
    JSparseB[902] = 1.0*dwdx426 + 1.0*dwdx427;
    JSparseB[903] = -1.0*dwdx20;
    JSparseB[904] = -1.0*dwdx394;
    JSparseB[905] = 1.0*dwdx428 + 1.0*dwdx429;
}