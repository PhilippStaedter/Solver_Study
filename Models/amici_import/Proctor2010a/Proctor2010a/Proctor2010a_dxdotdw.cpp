#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_Proctor2010a(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = 1.0;
    dxdotdw[1] = -1.0;
    dxdotdw[2] = 3.0;
    dxdotdw[3] = 3.0;
    dxdotdw[4] = 1.0;
    dxdotdw[5] = -1.0;
    dxdotdw[6] = 1.0;
    dxdotdw[7] = 1.0;
    dxdotdw[8] = -1.0;
    dxdotdw[9] = -1.0;
    dxdotdw[10] = -1.0;
    dxdotdw[11] = 1.0;
    dxdotdw[12] = 1.0;
    dxdotdw[13] = 1.0;
    dxdotdw[14] = -1.0;
    dxdotdw[15] = -1.0;
    dxdotdw[16] = -1.0;
    dxdotdw[17] = 1.0;
    dxdotdw[18] = 1.0;
    dxdotdw[19] = -1.0;
    dxdotdw[20] = -1.0;
    dxdotdw[21] = 1.0;
    dxdotdw[22] = 1.0;
    dxdotdw[23] = -1.0;
    dxdotdw[24] = -1.0;
    dxdotdw[25] = 1.0;
    dxdotdw[26] = 1.0;
    dxdotdw[27] = -1.0;
    dxdotdw[28] = -1.0;
    dxdotdw[29] = 1.0;
    dxdotdw[30] = 1.0;
    dxdotdw[31] = -1.0;
    dxdotdw[32] = -1.0;
    dxdotdw[33] = 1.0;
    dxdotdw[34] = 1.0;
    dxdotdw[35] = -1.0;
    dxdotdw[36] = -1.0;
    dxdotdw[37] = 1.0;
    dxdotdw[38] = 1.0;
    dxdotdw[39] = -1.0;
    dxdotdw[40] = -1.0;
    dxdotdw[41] = 1.0;
    dxdotdw[42] = 1.0;
    dxdotdw[43] = -1.0;
    dxdotdw[44] = -1.0;
    dxdotdw[45] = 1.0;
    dxdotdw[46] = 1.0;
    dxdotdw[47] = -1.0;
    dxdotdw[48] = -1.0;
    dxdotdw[49] = 1.0;
    dxdotdw[50] = 1.0;
    dxdotdw[51] = -1.0;
    dxdotdw[52] = -1.0;
    dxdotdw[53] = 1.0;
    dxdotdw[54] = -1.0;
    dxdotdw[55] = -1.0;
    dxdotdw[56] = 1.0;
    dxdotdw[57] = -1.0;
    dxdotdw[58] = -1.0;
    dxdotdw[59] = 1.0;
    dxdotdw[60] = -1.0;
    dxdotdw[61] = -1.0;
    dxdotdw[62] = 1.0;
    dxdotdw[63] = -1.0;
    dxdotdw[64] = -1.0;
    dxdotdw[65] = 1.0;
    dxdotdw[66] = -1.0;
    dxdotdw[67] = -1.0;
    dxdotdw[68] = 1.0;
    dxdotdw[69] = -1.0;
    dxdotdw[70] = -1.0;
    dxdotdw[71] = 1.0;
    dxdotdw[72] = -1.0;
    dxdotdw[73] = -1.0;
    dxdotdw[74] = 1.0;
    dxdotdw[75] = -1.0;
    dxdotdw[76] = -1.0;
    dxdotdw[77] = 1.0;
    dxdotdw[78] = 1.0;
    dxdotdw[79] = 1.0;
    dxdotdw[80] = -1.0;
    dxdotdw[81] = 1.0;
    dxdotdw[82] = 1.0;
    dxdotdw[83] = -1.0;
    dxdotdw[84] = 1.0;
    dxdotdw[85] = 1.0;
    dxdotdw[86] = -1.0;
    dxdotdw[87] = 1.0;
    dxdotdw[88] = 1.0;
    dxdotdw[89] = -1.0;
    dxdotdw[90] = 1.0;
    dxdotdw[91] = 1.0;
    dxdotdw[92] = -1.0;
    dxdotdw[93] = 1.0;
    dxdotdw[94] = 1.0;
    dxdotdw[95] = -1.0;
    dxdotdw[96] = 1.0;
    dxdotdw[97] = 1.0;
    dxdotdw[98] = -1.0;
    dxdotdw[99] = 1.0;
    dxdotdw[100] = 1.0;
    dxdotdw[101] = 1.0;
    dxdotdw[102] = -1.0;
    dxdotdw[103] = 1.0;
    dxdotdw[104] = -1.0;
    dxdotdw[105] = -1.0;
    dxdotdw[106] = 1.0;
    dxdotdw[107] = 1.0;
    dxdotdw[108] = -1.0;
    dxdotdw[109] = -1.0;
    dxdotdw[110] = 1.0;
    dxdotdw[111] = 1.0;
    dxdotdw[112] = -1.0;
    dxdotdw[113] = -1.0;
    dxdotdw[114] = 1.0;
    dxdotdw[115] = 1.0;
    dxdotdw[116] = -1.0;
    dxdotdw[117] = -1.0;
    dxdotdw[118] = 1.0;
    dxdotdw[119] = 1.0;
    dxdotdw[120] = -1.0;
    dxdotdw[121] = -1.0;
    dxdotdw[122] = 1.0;
    dxdotdw[123] = 1.0;
    dxdotdw[124] = 1.0;
    dxdotdw[125] = -1.0;
    dxdotdw[126] = 1.0;
    dxdotdw[127] = 1.0;
    dxdotdw[128] = -1.0;
    dxdotdw[129] = 1.0;
    dxdotdw[130] = 1.0;
    dxdotdw[131] = -1.0;
    dxdotdw[132] = 1.0;
    dxdotdw[133] = 1.0;
    dxdotdw[134] = -1.0;
    dxdotdw[135] = 1.0;
    dxdotdw[136] = 4.0;
    dxdotdw[137] = 1.0;
    dxdotdw[138] = -1.0;
    dxdotdw[139] = 4.0;
    dxdotdw[140] = 1.0;
    dxdotdw[141] = -1.0;
    dxdotdw[142] = 5.0;
    dxdotdw[143] = 1.0;
    dxdotdw[144] = -1.0;
    dxdotdw[145] = 6.0;
    dxdotdw[146] = 1.0;
    dxdotdw[147] = -1.0;
    dxdotdw[148] = 7.0;
    dxdotdw[149] = 1.0;
    dxdotdw[150] = -1.0;
    dxdotdw[151] = 8.0;
    dxdotdw[152] = 1.0;
    dxdotdw[153] = -1.0;
    dxdotdw[154] = -2.0;
    dxdotdw[155] = 1.0;
    dxdotdw[156] = -1.0;
    dxdotdw[157] = -1.0;
    dxdotdw[158] = 1.0;
    dxdotdw[159] = -1.0;
    dxdotdw[160] = -1.0;
    dxdotdw[161] = 1.0;
    dxdotdw[162] = -1.0;
    dxdotdw[163] = -1.0;
    dxdotdw[164] = 1.0;
    dxdotdw[165] = -1.0;
    dxdotdw[166] = -1.0;
    dxdotdw[167] = 1.0;
    dxdotdw[168] = 1.0;
    dxdotdw[169] = 1.0;
    dxdotdw[170] = -1.0;
    dxdotdw[171] = 1.0;
    dxdotdw[172] = 1.0;
    dxdotdw[173] = -1.0;
    dxdotdw[174] = 1.0;
    dxdotdw[175] = 1.0;
    dxdotdw[176] = -1.0;
    dxdotdw[177] = 1.0;
    dxdotdw[178] = 1.0;
    dxdotdw[179] = -1.0;
    dxdotdw[180] = 2.0;
    dxdotdw[181] = -1.0;
    dxdotdw[182] = -1.0;
    dxdotdw[183] = -1.0;
    dxdotdw[184] = 1.0;
    dxdotdw[185] = 7.0;
    dxdotdw[186] = -1.0;
    dxdotdw[187] = 1.0;
    dxdotdw[188] = 1.0;
    dxdotdw[189] = -1.0;
    dxdotdw[190] = 1.0;
    dxdotdw[191] = 1.0;
    dxdotdw[192] = 1.0;
    dxdotdw[193] = -1.0;
    dxdotdw[194] = 1.0;
    dxdotdw[195] = 1.0;
    dxdotdw[196] = 1.0;
    dxdotdw[197] = 1.0;
    dxdotdw[198] = -1.0;
    dxdotdw[199] = 1.0;
    dxdotdw[200] = 2.0;
    dxdotdw[201] = 1.0;
    dxdotdw[202] = 1.0;
    dxdotdw[203] = -1.0;
    dxdotdw[204] = 1.0;
    dxdotdw[205] = 3.0;
    dxdotdw[206] = 1.0;
    dxdotdw[207] = 1.0;
    dxdotdw[208] = -1.0;
    dxdotdw[209] = 1.0;
    dxdotdw[210] = 4.0;
    dxdotdw[211] = 1.0;
    dxdotdw[212] = 1.0;
    dxdotdw[213] = -1.0;
    dxdotdw[214] = 1.0;
    dxdotdw[215] = 5.0;
    dxdotdw[216] = 1.0;
    dxdotdw[217] = 1.0;
    dxdotdw[218] = -1.0;
    dxdotdw[219] = 1.0;
    dxdotdw[220] = 6.0;
    dxdotdw[221] = 1.0;
    dxdotdw[222] = 1.0;
    dxdotdw[223] = -1.0;
    dxdotdw[224] = 1.0;
    dxdotdw[225] = 7.0;
    dxdotdw[226] = 1.0;
    dxdotdw[227] = 1.0;
    dxdotdw[228] = -1.0;
    dxdotdw[229] = 1.0;
    dxdotdw[230] = 8.0;
    dxdotdw[231] = 1.0;
    dxdotdw[232] = -1.0;
    dxdotdw[233] = 1.0;
    dxdotdw[234] = 1.0;
    dxdotdw[235] = 1.0;
    dxdotdw[236] = 1.0;
    dxdotdw[237] = 1.0;
    dxdotdw[238] = -1.0;
    dxdotdw[239] = 1.0;
    dxdotdw[240] = 2.0;
    dxdotdw[241] = 1.0;
    dxdotdw[242] = 1.0;
    dxdotdw[243] = 1.0;
    dxdotdw[244] = -1.0;
    dxdotdw[245] = 1.0;
    dxdotdw[246] = 3.0;
    dxdotdw[247] = 1.0;
    dxdotdw[248] = 1.0;
    dxdotdw[249] = 1.0;
    dxdotdw[250] = -1.0;
    dxdotdw[251] = 1.0;
    dxdotdw[252] = 4.0;
    dxdotdw[253] = 1.0;
    dxdotdw[254] = 1.0;
    dxdotdw[255] = 1.0;
    dxdotdw[256] = -1.0;
    dxdotdw[257] = 1.0;
    dxdotdw[258] = 5.0;
    dxdotdw[259] = 1.0;
    dxdotdw[260] = 1.0;
    dxdotdw[261] = 1.0;
    dxdotdw[262] = -1.0;
    dxdotdw[263] = 1.0;
    dxdotdw[264] = 6.0;
    dxdotdw[265] = 1.0;
    dxdotdw[266] = 1.0;
    dxdotdw[267] = 1.0;
    dxdotdw[268] = -1.0;
    dxdotdw[269] = 1.0;
    dxdotdw[270] = 7.0;
    dxdotdw[271] = 1.0;
    dxdotdw[272] = 1.0;
    dxdotdw[273] = 1.0;
    dxdotdw[274] = -1.0;
    dxdotdw[275] = 1.0;
    dxdotdw[276] = 8.0;
    dxdotdw[277] = 1.0;
    dxdotdw[278] = 1.0;
    dxdotdw[279] = 1.0;
    dxdotdw[280] = -1.0;
    dxdotdw[281] = -1.0;
    dxdotdw[282] = 1.0;
    dxdotdw[283] = -1.0;
    dxdotdw[284] = -1.0;
    dxdotdw[285] = 1.0;
    dxdotdw[286] = -1.0;
    dxdotdw[287] = -1.0;
    dxdotdw[288] = 1.0;
    dxdotdw[289] = -1.0;
    dxdotdw[290] = -1.0;
    dxdotdw[291] = 1.0;
    dxdotdw[292] = -1.0;
    dxdotdw[293] = -1.0;
    dxdotdw[294] = 1.0;
    dxdotdw[295] = 1.0;
    dxdotdw[296] = 1.0;
    dxdotdw[297] = 1.0;
    dxdotdw[298] = 1.0;
    dxdotdw[299] = 1.0;
    dxdotdw[300] = 1.0;
    dxdotdw[301] = -1.0;
    dxdotdw[302] = -1.0;
    dxdotdw[303] = 1.0;
    dxdotdw[304] = 1.0;
    dxdotdw[305] = -1.0;
    dxdotdw[306] = -1.0;
    dxdotdw[307] = -1.0;
    dxdotdw[308] = 1.0;
    dxdotdw[309] = -1.0;
    dxdotdw[310] = 1.0;
    dxdotdw[311] = -1.0;
    dxdotdw[312] = 1.0;
    dxdotdw[313] = -1.0;
    dxdotdw[314] = -1.0;
    dxdotdw[315] = -1.0;
    dxdotdw[316] = -1.0;
    dxdotdw[317] = 1.0;
    dxdotdw[318] = 1.0;
    dxdotdw[319] = 1.0;
    dxdotdw[320] = -1.0;
    dxdotdw[321] = -1.0;
    dxdotdw[322] = -1.0;
    dxdotdw[323] = 1.0;
    dxdotdw[324] = 1.0;
    dxdotdw[325] = 1.0;
    dxdotdw[326] = -1.0;
    dxdotdw[327] = 1.0;
    dxdotdw[328] = -1.0;
    dxdotdw[329] = 1.0;
    dxdotdw[330] = 1.0;
    dxdotdw[331] = -1.0;
    dxdotdw[332] = -1.0;
    dxdotdw[333] = -1.0;
    dxdotdw[334] = 1.0;
    dxdotdw[335] = 1.0;
    dxdotdw[336] = 1.0;
    dxdotdw[337] = -1.0;
    dxdotdw[338] = 1.0;
    dxdotdw[339] = -1.0;
    dxdotdw[340] = -1.0;
    dxdotdw[341] = 1.0;
    dxdotdw[342] = 1.0;
    dxdotdw[343] = -1.0;
    dxdotdw[344] = -1.0;
    dxdotdw[345] = 1.0;
    dxdotdw[346] = 1.0;
    dxdotdw[347] = -1.0;
    dxdotdw[348] = -1.0;
    dxdotdw[349] = 1.0;
    dxdotdw[350] = 1.0;
    dxdotdw[351] = -1.0;
    dxdotdw[352] = -1.0;
    dxdotdw[353] = 1.0;
    dxdotdw[354] = 1.0;
    dxdotdw[355] = -1.0;
    dxdotdw[356] = -1.0;
    dxdotdw[357] = 1.0;
    dxdotdw[358] = 1.0;
    dxdotdw[359] = -1.0;
    dxdotdw[360] = -1.0;
    dxdotdw[361] = 1.0;
    dxdotdw[362] = 1.0;
    dxdotdw[363] = -1.0;
    dxdotdw[364] = -1.0;
    dxdotdw[365] = 1.0;
    dxdotdw[366] = 1.0;
    dxdotdw[367] = -1.0;
    dxdotdw[368] = -1.0;
    dxdotdw[369] = 1.0;
    dxdotdw[370] = -1.0;
    dxdotdw[371] = -1.0;
    dxdotdw[372] = 1.0;
    dxdotdw[373] = -1.0;
    dxdotdw[374] = -1.0;
    dxdotdw[375] = 1.0;
    dxdotdw[376] = -1.0;
    dxdotdw[377] = -1.0;
    dxdotdw[378] = 1.0;
    dxdotdw[379] = -1.0;
    dxdotdw[380] = -1.0;
    dxdotdw[381] = 1.0;
    dxdotdw[382] = -1.0;
    dxdotdw[383] = -1.0;
    dxdotdw[384] = 1.0;
    dxdotdw[385] = -1.0;
    dxdotdw[386] = -1.0;
    dxdotdw[387] = 1.0;
    dxdotdw[388] = -1.0;
    dxdotdw[389] = -1.0;
    dxdotdw[390] = 1.0;
    dxdotdw[391] = -1.0;
    dxdotdw[392] = -1.0;
    dxdotdw[393] = 1.0;
    dxdotdw[394] = 1.0;
    dxdotdw[395] = 1.0;
    dxdotdw[396] = -1.0;
    dxdotdw[397] = 1.0;
    dxdotdw[398] = 1.0;
    dxdotdw[399] = -1.0;
    dxdotdw[400] = 1.0;
    dxdotdw[401] = 1.0;
    dxdotdw[402] = -1.0;
    dxdotdw[403] = 1.0;
    dxdotdw[404] = 1.0;
    dxdotdw[405] = -1.0;
    dxdotdw[406] = 1.0;
    dxdotdw[407] = 1.0;
    dxdotdw[408] = -1.0;
    dxdotdw[409] = 1.0;
    dxdotdw[410] = 1.0;
    dxdotdw[411] = -1.0;
    dxdotdw[412] = 1.0;
    dxdotdw[413] = 1.0;
    dxdotdw[414] = -1.0;
    dxdotdw[415] = 1.0;
    dxdotdw[416] = 1.0;
    dxdotdw[417] = 1.0;
    dxdotdw[418] = -1.0;
    dxdotdw[419] = -1.0;
    dxdotdw[420] = 1.0;
    dxdotdw[421] = -1.0;
    dxdotdw[422] = 1.0;
    dxdotdw[423] = -1.0;
    dxdotdw[424] = 1.0;
    dxdotdw[425] = -1.0;
    dxdotdw[426] = 1.0;
    dxdotdw[427] = -1.0;
    dxdotdw[428] = 1.0;
    dxdotdw[429] = -1.0;
    dxdotdw[430] = 1.0;
    dxdotdw[431] = -1.0;
    dxdotdw[432] = 1.0;
    dxdotdw[433] = -1.0;
    dxdotdw[434] = 1.0;
    dxdotdw[435] = -1.0;
    dxdotdw[436] = 1.0;
    dxdotdw[437] = -1.0;
    dxdotdw[438] = 1.0;
    dxdotdw[439] = 1.0;
    dxdotdw[440] = 1.0;
    dxdotdw[441] = -1.0;
    dxdotdw[442] = 1.0;
    dxdotdw[443] = 1.0;
    dxdotdw[444] = -1.0;
    dxdotdw[445] = 1.0;
    dxdotdw[446] = 1.0;
    dxdotdw[447] = -1.0;
    dxdotdw[448] = 1.0;
    dxdotdw[449] = 1.0;
    dxdotdw[450] = -1.0;
    dxdotdw[451] = 4.0;
    dxdotdw[452] = 1.0;
    dxdotdw[453] = 1.0;
    dxdotdw[454] = -1.0;
    dxdotdw[455] = 4.0;
    dxdotdw[456] = 1.0;
    dxdotdw[457] = -1.0;
    dxdotdw[458] = 5.0;
    dxdotdw[459] = 1.0;
    dxdotdw[460] = -1.0;
    dxdotdw[461] = 6.0;
    dxdotdw[462] = 1.0;
    dxdotdw[463] = -1.0;
    dxdotdw[464] = 7.0;
    dxdotdw[465] = 1.0;
    dxdotdw[466] = -1.0;
    dxdotdw[467] = 8.0;
    dxdotdw[468] = 1.0;
    dxdotdw[469] = -1.0;
    dxdotdw[470] = 1.0;
    dxdotdw[471] = -1.0;
    dxdotdw[472] = -1.0;
    dxdotdw[473] = 1.0;
    dxdotdw[474] = 1.0;
    dxdotdw[475] = -1.0;
    dxdotdw[476] = -1.0;
    dxdotdw[477] = -1.0;
    dxdotdw[478] = 1.0;
    dxdotdw[479] = 1.0;
    dxdotdw[480] = -1.0;
    dxdotdw[481] = -1.0;
    dxdotdw[482] = 1.0;
    dxdotdw[483] = -1.0;
    dxdotdw[484] = -1.0;
    dxdotdw[485] = 1.0;
    dxdotdw[486] = 1.0;
    dxdotdw[487] = 1.0;
    dxdotdw[488] = -1.0;
    dxdotdw[489] = 1.0;
    dxdotdw[490] = -1.0;
    dxdotdw[491] = -1.0;
    dxdotdw[492] = 1.0;
    dxdotdw[493] = 1.0;
    dxdotdw[494] = -1.0;
    dxdotdw[495] = -1.0;
    dxdotdw[496] = 1.0;
    dxdotdw[497] = 1.0;
    dxdotdw[498] = -1.0;
    dxdotdw[499] = -1.0;
    dxdotdw[500] = 1.0;
    dxdotdw[501] = 1.0;
    dxdotdw[502] = -1.0;
    dxdotdw[503] = -1.0;
    dxdotdw[504] = 1.0;
    dxdotdw[505] = 1.0;
    dxdotdw[506] = -1.0;
    dxdotdw[507] = -1.0;
    dxdotdw[508] = 1.0;
    dxdotdw[509] = 1.0;
    dxdotdw[510] = -1.0;
    dxdotdw[511] = -1.0;
    dxdotdw[512] = 1.0;
    dxdotdw[513] = 1.0;
    dxdotdw[514] = -1.0;
    dxdotdw[515] = -1.0;
    dxdotdw[516] = 1.0;
    dxdotdw[517] = 1.0;
    dxdotdw[518] = -1.0;
    dxdotdw[519] = -1.0;
    dxdotdw[520] = 1.0;
    dxdotdw[521] = -1.0;
    dxdotdw[522] = -1.0;
    dxdotdw[523] = 1.0;
    dxdotdw[524] = -1.0;
    dxdotdw[525] = -1.0;
    dxdotdw[526] = 1.0;
    dxdotdw[527] = -1.0;
    dxdotdw[528] = -1.0;
    dxdotdw[529] = 1.0;
    dxdotdw[530] = -1.0;
    dxdotdw[531] = -1.0;
    dxdotdw[532] = 1.0;
    dxdotdw[533] = -1.0;
    dxdotdw[534] = -1.0;
    dxdotdw[535] = 1.0;
    dxdotdw[536] = -1.0;
    dxdotdw[537] = -1.0;
    dxdotdw[538] = 1.0;
    dxdotdw[539] = -1.0;
    dxdotdw[540] = -1.0;
    dxdotdw[541] = 1.0;
    dxdotdw[542] = -1.0;
    dxdotdw[543] = -1.0;
    dxdotdw[544] = 1.0;
    dxdotdw[545] = 1.0;
    dxdotdw[546] = 1.0;
    dxdotdw[547] = -1.0;
    dxdotdw[548] = 1.0;
    dxdotdw[549] = 1.0;
    dxdotdw[550] = -1.0;
    dxdotdw[551] = 1.0;
    dxdotdw[552] = 1.0;
    dxdotdw[553] = -1.0;
    dxdotdw[554] = 1.0;
    dxdotdw[555] = 1.0;
    dxdotdw[556] = -1.0;
    dxdotdw[557] = 1.0;
    dxdotdw[558] = 1.0;
    dxdotdw[559] = -1.0;
    dxdotdw[560] = 1.0;
    dxdotdw[561] = 1.0;
    dxdotdw[562] = -1.0;
    dxdotdw[563] = 1.0;
    dxdotdw[564] = 1.0;
    dxdotdw[565] = -1.0;
    dxdotdw[566] = 1.0;
    dxdotdw[567] = 1.0;
    dxdotdw[568] = 1.0;
    dxdotdw[569] = -1.0;
    dxdotdw[570] = -1.0;
    dxdotdw[571] = 1.0;
    dxdotdw[572] = -1.0;
    dxdotdw[573] = 1.0;
    dxdotdw[574] = -1.0;
    dxdotdw[575] = 1.0;
    dxdotdw[576] = -1.0;
    dxdotdw[577] = 1.0;
    dxdotdw[578] = -1.0;
    dxdotdw[579] = 1.0;
    dxdotdw[580] = -1.0;
    dxdotdw[581] = 1.0;
    dxdotdw[582] = -1.0;
    dxdotdw[583] = 1.0;
    dxdotdw[584] = -1.0;
    dxdotdw[585] = 1.0;
    dxdotdw[586] = -1.0;
    dxdotdw[587] = 1.0;
    dxdotdw[588] = -1.0;
    dxdotdw[589] = 1.0;
    dxdotdw[590] = 1.0;
    dxdotdw[591] = 1.0;
    dxdotdw[592] = -1.0;
    dxdotdw[593] = 1.0;
    dxdotdw[594] = 1.0;
    dxdotdw[595] = -1.0;
    dxdotdw[596] = 1.0;
    dxdotdw[597] = 1.0;
    dxdotdw[598] = -1.0;
    dxdotdw[599] = 1.0;
    dxdotdw[600] = 1.0;
    dxdotdw[601] = -1.0;
    dxdotdw[602] = 4.0;
    dxdotdw[603] = 1.0;
    dxdotdw[604] = 1.0;
    dxdotdw[605] = -1.0;
    dxdotdw[606] = 4.0;
    dxdotdw[607] = 1.0;
    dxdotdw[608] = -1.0;
    dxdotdw[609] = 5.0;
    dxdotdw[610] = 1.0;
    dxdotdw[611] = -1.0;
    dxdotdw[612] = 6.0;
    dxdotdw[613] = 1.0;
    dxdotdw[614] = -1.0;
    dxdotdw[615] = 7.0;
    dxdotdw[616] = 1.0;
    dxdotdw[617] = -1.0;
    dxdotdw[618] = 8.0;
    dxdotdw[619] = 1.0;
    dxdotdw[620] = -1.0;
    dxdotdw[621] = -2.0;
    dxdotdw[622] = 1.0;
    dxdotdw[623] = -1.0;
    dxdotdw[624] = -1.0;
    dxdotdw[625] = 1.0;
    dxdotdw[626] = -1.0;
    dxdotdw[627] = -1.0;
    dxdotdw[628] = 1.0;
    dxdotdw[629] = -1.0;
    dxdotdw[630] = -1.0;
    dxdotdw[631] = 1.0;
    dxdotdw[632] = -1.0;
    dxdotdw[633] = -1.0;
    dxdotdw[634] = 1.0;
    dxdotdw[635] = 1.0;
    dxdotdw[636] = 1.0;
    dxdotdw[637] = -1.0;
    dxdotdw[638] = 1.0;
    dxdotdw[639] = 1.0;
    dxdotdw[640] = -1.0;
    dxdotdw[641] = 1.0;
    dxdotdw[642] = 1.0;
    dxdotdw[643] = -1.0;
    dxdotdw[644] = 1.0;
    dxdotdw[645] = 1.0;
    dxdotdw[646] = -1.0;
    dxdotdw[647] = 2.0;
    dxdotdw[648] = -1.0;
    dxdotdw[649] = -1.0;
    dxdotdw[650] = 1.0;
    dxdotdw[651] = -1.0;
    dxdotdw[652] = -1.0;
    dxdotdw[653] = 1.0;
    dxdotdw[654] = -1.0;
    dxdotdw[655] = -1.0;
    dxdotdw[656] = 1.0;
    dxdotdw[657] = -1.0;
    dxdotdw[658] = -1.0;
    dxdotdw[659] = 1.0;
    dxdotdw[660] = -1.0;
    dxdotdw[661] = -1.0;
    dxdotdw[662] = 1.0;
    dxdotdw[663] = -1.0;
    dxdotdw[664] = 1.0;
    dxdotdw[665] = -1.0;
    dxdotdw[666] = -1.0;
    dxdotdw[667] = 7.0;
    dxdotdw[668] = 1.0;
    dxdotdw[669] = -1.0;
    dxdotdw[670] = 1.0;
    dxdotdw[671] = 1.0;
    dxdotdw[672] = 1.0;
    dxdotdw[673] = 1.0;
    dxdotdw[674] = 1.0;
    dxdotdw[675] = 1.0;
    dxdotdw[676] = -2.0;
    dxdotdw[677] = 1.0;
    dxdotdw[678] = -1.0;
    dxdotdw[679] = -1.0;
    dxdotdw[680] = 1.0;
    dxdotdw[681] = -1.0;
    dxdotdw[682] = -1.0;
    dxdotdw[683] = 1.0;
    dxdotdw[684] = -1.0;
    dxdotdw[685] = -1.0;
    dxdotdw[686] = 1.0;
    dxdotdw[687] = -1.0;
    dxdotdw[688] = -1.0;
    dxdotdw[689] = 1.0;
    dxdotdw[690] = 1.0;
    dxdotdw[691] = 1.0;
    dxdotdw[692] = -1.0;
    dxdotdw[693] = 1.0;
    dxdotdw[694] = 1.0;
    dxdotdw[695] = -1.0;
    dxdotdw[696] = 1.0;
    dxdotdw[697] = 1.0;
    dxdotdw[698] = -1.0;
    dxdotdw[699] = 1.0;
    dxdotdw[700] = 1.0;
    dxdotdw[701] = -1.0;
    dxdotdw[702] = 2.0;
    dxdotdw[703] = -1.0;
    dxdotdw[704] = -1.0;
    dxdotdw[705] = 1.0;
    dxdotdw[706] = -1.0;
    dxdotdw[707] = -1.0;
    dxdotdw[708] = 1.0;
    dxdotdw[709] = -1.0;
    dxdotdw[710] = -1.0;
    dxdotdw[711] = 1.0;
    dxdotdw[712] = -1.0;
    dxdotdw[713] = -1.0;
    dxdotdw[714] = 1.0;
    dxdotdw[715] = -1.0;
    dxdotdw[716] = -1.0;
    dxdotdw[717] = 1.0;
    dxdotdw[718] = -1.0;
    dxdotdw[719] = 1.0;
    dxdotdw[720] = -1.0;
    dxdotdw[721] = -1.0;
    dxdotdw[722] = 7.0;
    dxdotdw[723] = 1.0;
    dxdotdw[724] = -1.0;
    dxdotdw[725] = 1.0;
    dxdotdw[726] = 1.0;
    dxdotdw[727] = -1.0;
    dxdotdw[728] = 1.0;
    dxdotdw[729] = 1.0;
    dxdotdw[730] = 1.0;
    dxdotdw[731] = -1.0;
    dxdotdw[732] = 1.0;
    dxdotdw[733] = 1.0;
    dxdotdw[734] = 1.0;
    dxdotdw[735] = 1.0;
    dxdotdw[736] = -1.0;
    dxdotdw[737] = 1.0;
    dxdotdw[738] = 1.0;
    dxdotdw[739] = 2.0;
    dxdotdw[740] = 1.0;
    dxdotdw[741] = -1.0;
    dxdotdw[742] = 1.0;
    dxdotdw[743] = 1.0;
    dxdotdw[744] = 3.0;
    dxdotdw[745] = 1.0;
    dxdotdw[746] = -1.0;
    dxdotdw[747] = 1.0;
    dxdotdw[748] = 1.0;
    dxdotdw[749] = 4.0;
    dxdotdw[750] = 1.0;
    dxdotdw[751] = -1.0;
    dxdotdw[752] = 1.0;
    dxdotdw[753] = 1.0;
    dxdotdw[754] = 5.0;
    dxdotdw[755] = 1.0;
    dxdotdw[756] = -1.0;
    dxdotdw[757] = 1.0;
    dxdotdw[758] = 1.0;
    dxdotdw[759] = 6.0;
    dxdotdw[760] = 1.0;
    dxdotdw[761] = -1.0;
    dxdotdw[762] = 1.0;
    dxdotdw[763] = 1.0;
    dxdotdw[764] = 7.0;
    dxdotdw[765] = 1.0;
    dxdotdw[766] = -1.0;
    dxdotdw[767] = 1.0;
    dxdotdw[768] = 1.0;
    dxdotdw[769] = 8.0;
    dxdotdw[770] = 1.0;
    dxdotdw[771] = -1.0;
    dxdotdw[772] = 1.0;
    dxdotdw[773] = 1.0;
    dxdotdw[774] = 1.0;
    dxdotdw[775] = 1.0;
    dxdotdw[776] = -1.0;
    dxdotdw[777] = 1.0;
    dxdotdw[778] = 1.0;
    dxdotdw[779] = 2.0;
    dxdotdw[780] = 1.0;
    dxdotdw[781] = -1.0;
    dxdotdw[782] = 1.0;
    dxdotdw[783] = 1.0;
    dxdotdw[784] = 3.0;
    dxdotdw[785] = 1.0;
    dxdotdw[786] = -1.0;
    dxdotdw[787] = 1.0;
    dxdotdw[788] = 1.0;
    dxdotdw[789] = 4.0;
    dxdotdw[790] = 1.0;
    dxdotdw[791] = -1.0;
    dxdotdw[792] = 1.0;
    dxdotdw[793] = 1.0;
    dxdotdw[794] = 5.0;
    dxdotdw[795] = 1.0;
    dxdotdw[796] = -1.0;
    dxdotdw[797] = 1.0;
    dxdotdw[798] = 1.0;
    dxdotdw[799] = 6.0;
    dxdotdw[800] = 1.0;
    dxdotdw[801] = -1.0;
    dxdotdw[802] = 1.0;
    dxdotdw[803] = 1.0;
    dxdotdw[804] = 7.0;
    dxdotdw[805] = 1.0;
    dxdotdw[806] = -1.0;
    dxdotdw[807] = 1.0;
    dxdotdw[808] = 1.0;
    dxdotdw[809] = 8.0;
    dxdotdw[810] = 1.0;
    dxdotdw[811] = 1.0;
    dxdotdw[812] = 1.0;
    dxdotdw[813] = 1.0;
    dxdotdw[814] = 1.0;
    dxdotdw[815] = -2.0;
    dxdotdw[816] = 1.0;
    dxdotdw[817] = -1.0;
    dxdotdw[818] = -1.0;
    dxdotdw[819] = 1.0;
    dxdotdw[820] = -1.0;
    dxdotdw[821] = -1.0;
    dxdotdw[822] = 1.0;
    dxdotdw[823] = -1.0;
    dxdotdw[824] = -1.0;
    dxdotdw[825] = 1.0;
    dxdotdw[826] = -1.0;
    dxdotdw[827] = -1.0;
    dxdotdw[828] = 1.0;
    dxdotdw[829] = 1.0;
    dxdotdw[830] = 1.0;
    dxdotdw[831] = -1.0;
    dxdotdw[832] = 1.0;
    dxdotdw[833] = 1.0;
    dxdotdw[834] = -1.0;
    dxdotdw[835] = 1.0;
    dxdotdw[836] = 1.0;
    dxdotdw[837] = -1.0;
    dxdotdw[838] = 1.0;
    dxdotdw[839] = 1.0;
    dxdotdw[840] = -1.0;
    dxdotdw[841] = 2.0;
    dxdotdw[842] = -1.0;
    dxdotdw[843] = -1.0;
    dxdotdw[844] = 1.0;
    dxdotdw[845] = -1.0;
    dxdotdw[846] = -1.0;
    dxdotdw[847] = 1.0;
    dxdotdw[848] = -1.0;
    dxdotdw[849] = -1.0;
    dxdotdw[850] = 1.0;
    dxdotdw[851] = -1.0;
    dxdotdw[852] = -1.0;
    dxdotdw[853] = 1.0;
    dxdotdw[854] = -1.0;
    dxdotdw[855] = -1.0;
    dxdotdw[856] = 1.0;
    dxdotdw[857] = -1.0;
    dxdotdw[858] = 1.0;
    dxdotdw[859] = -1.0;
    dxdotdw[860] = -1.0;
    dxdotdw[861] = 7.0;
    dxdotdw[862] = 1.0;
    dxdotdw[863] = -1.0;
    dxdotdw[864] = 1.0;
    dxdotdw[865] = 1.0;
    dxdotdw[866] = -1.0;
    dxdotdw[867] = 1.0;
    dxdotdw[868] = 1.0;
    dxdotdw[869] = 1.0;
    dxdotdw[870] = 1.0;
    dxdotdw[871] = 1.0;
    dxdotdw[872] = 1.0;
    dxdotdw[873] = -2.0;
    dxdotdw[874] = 1.0;
    dxdotdw[875] = -1.0;
    dxdotdw[876] = -1.0;
    dxdotdw[877] = 1.0;
    dxdotdw[878] = -1.0;
    dxdotdw[879] = -1.0;
    dxdotdw[880] = 1.0;
    dxdotdw[881] = -1.0;
    dxdotdw[882] = -1.0;
    dxdotdw[883] = 1.0;
    dxdotdw[884] = -1.0;
    dxdotdw[885] = -1.0;
    dxdotdw[886] = 1.0;
    dxdotdw[887] = 1.0;
    dxdotdw[888] = 1.0;
    dxdotdw[889] = -1.0;
    dxdotdw[890] = 1.0;
    dxdotdw[891] = 1.0;
    dxdotdw[892] = -1.0;
    dxdotdw[893] = 1.0;
    dxdotdw[894] = 1.0;
    dxdotdw[895] = -1.0;
    dxdotdw[896] = 1.0;
    dxdotdw[897] = 1.0;
    dxdotdw[898] = -1.0;
    dxdotdw[899] = 2.0;
    dxdotdw[900] = -1.0;
    dxdotdw[901] = -1.0;
    dxdotdw[902] = 1.0;
    dxdotdw[903] = -1.0;
    dxdotdw[904] = -1.0;
    dxdotdw[905] = 1.0;
    dxdotdw[906] = -1.0;
    dxdotdw[907] = -1.0;
    dxdotdw[908] = 1.0;
    dxdotdw[909] = -1.0;
    dxdotdw[910] = -1.0;
    dxdotdw[911] = 1.0;
    dxdotdw[912] = -1.0;
    dxdotdw[913] = -1.0;
    dxdotdw[914] = 1.0;
    dxdotdw[915] = -1.0;
    dxdotdw[916] = 1.0;
    dxdotdw[917] = -1.0;
    dxdotdw[918] = -1.0;
    dxdotdw[919] = 7.0;
    dxdotdw[920] = 1.0;
    dxdotdw[921] = -1.0;
    dxdotdw[922] = 1.0;
    dxdotdw[923] = 1.0;
    dxdotdw[924] = -1.0;
    dxdotdw[925] = 1.0;
    dxdotdw[926] = 1.0;
    dxdotdw[927] = -1.0;
    dxdotdw[928] = 1.0;
    dxdotdw[929] = 1.0;
    dxdotdw[930] = 1.0;
    dxdotdw[931] = -1.0;
    dxdotdw[932] = 2.0;
    dxdotdw[933] = 1.0;
    dxdotdw[934] = 1.0;
    dxdotdw[935] = -1.0;
    dxdotdw[936] = 3.0;
    dxdotdw[937] = 1.0;
    dxdotdw[938] = 1.0;
    dxdotdw[939] = -1.0;
    dxdotdw[940] = 4.0;
    dxdotdw[941] = 1.0;
    dxdotdw[942] = 1.0;
    dxdotdw[943] = -1.0;
    dxdotdw[944] = 5.0;
    dxdotdw[945] = 1.0;
    dxdotdw[946] = 1.0;
    dxdotdw[947] = -1.0;
    dxdotdw[948] = 6.0;
    dxdotdw[949] = 1.0;
    dxdotdw[950] = 1.0;
    dxdotdw[951] = -1.0;
    dxdotdw[952] = 7.0;
    dxdotdw[953] = 1.0;
    dxdotdw[954] = 1.0;
    dxdotdw[955] = -1.0;
    dxdotdw[956] = 8.0;
    dxdotdw[957] = 1.0;
    dxdotdw[958] = 1.0;
    dxdotdw[959] = -1.0;
    dxdotdw[960] = 1.0;
    dxdotdw[961] = 1.0;
    dxdotdw[962] = 1.0;
    dxdotdw[963] = 1.0;
    dxdotdw[964] = -1.0;
    dxdotdw[965] = 2.0;
    dxdotdw[966] = 1.0;
    dxdotdw[967] = 1.0;
    dxdotdw[968] = 1.0;
    dxdotdw[969] = -1.0;
    dxdotdw[970] = 3.0;
    dxdotdw[971] = 1.0;
    dxdotdw[972] = 1.0;
    dxdotdw[973] = 1.0;
    dxdotdw[974] = -1.0;
    dxdotdw[975] = 4.0;
    dxdotdw[976] = 1.0;
    dxdotdw[977] = 1.0;
    dxdotdw[978] = 1.0;
    dxdotdw[979] = -1.0;
    dxdotdw[980] = 5.0;
    dxdotdw[981] = 1.0;
    dxdotdw[982] = 1.0;
    dxdotdw[983] = 1.0;
    dxdotdw[984] = -1.0;
    dxdotdw[985] = 6.0;
    dxdotdw[986] = 1.0;
    dxdotdw[987] = 1.0;
    dxdotdw[988] = 1.0;
    dxdotdw[989] = -1.0;
    dxdotdw[990] = 7.0;
    dxdotdw[991] = 1.0;
    dxdotdw[992] = 1.0;
    dxdotdw[993] = 1.0;
    dxdotdw[994] = -1.0;
    dxdotdw[995] = 8.0;
    dxdotdw[996] = 1.0;
    dxdotdw[997] = 1.0;
    dxdotdw[998] = 1.0;
    dxdotdw[999] = 1.0;
    dxdotdw[1000] = 1.0;
    dxdotdw[1001] = 1.0;
    dxdotdw[1002] = 1.0;
    dxdotdw[1003] = 1.0;
    dxdotdw[1004] = -1.0;
}