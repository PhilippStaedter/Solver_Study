#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_Pritchard2002(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[1] = -1.0*dwdx0;
    JB[26] = -1.0*dwdx1 + 1.0*dwdx2;
    JB[27] = 1.0*dwdx2;
    JB[28] = -1.0*dwdx2;
    JB[29] = -1.0*dwdx2;
    JB[51] = 1.0*dwdx3;
    JB[52] = 1.0*dwdx3 + 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7 - 1.0*dwdx8;
    JB[53] = -1.0*dwdx3;
    JB[54] = -1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7 + 2.0*dwdx8;
    JB[55] = 1.0*dwdx4;
    JB[56] = -1.0*dwdx4;
    JB[57] = -1.0*dwdx8;
    JB[62] = 1.0*dwdx5;
    JB[64] = -1.0*dwdx5;
    JB[66] = 1.0*dwdx6;
    JB[67] = -1.0*dwdx6;
    JB[76] = 1.0*dwdx9;
    JB[77] = 1.0*dwdx9;
    JB[78] = 1.0*dwdx10 - 1.0*dwdx9;
    JB[79] = -1.0*dwdx9;
    JB[80] = -1.0*dwdx10;
    JB[101] = 1.0*dwdx11;
    JB[102] = 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14;
    JB[103] = -1.0*dwdx11;
    JB[104] = -1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 2.0*dwdx14;
    JB[107] = -1.0*dwdx14;
    JB[112] = 1.0*dwdx12;
    JB[114] = -1.0*dwdx12;
    JB[116] = 1.0*dwdx13;
    JB[117] = -1.0*dwdx13;
    JB[127] = 1.0*dwdx16;
    JB[128] = 1.0*dwdx15;
    JB[129] = -1.0*dwdx16;
    JB[130] = -1.0*dwdx15 + 1.0*dwdx16;
    JB[131] = -1.0*dwdx16;
    JB[152] = 1.0*dwdx17;
    JB[154] = -1.0*dwdx17;
    JB[155] = 1.0*dwdx17;
    JB[156] = -1.0*dwdx17 + 1.0*dwdx18;
    JB[159] = -1.0*dwdx18;
    JB[160] = -1.0*dwdx18;
    JB[177] = 1.0*dwdx19 - 1.0*dwdx20;
    JB[179] = -1.0*dwdx19 + 2.0*dwdx20;
    JB[180] = 1.0*dwdx19;
    JB[181] = -1.0*dwdx19;
    JB[182] = -1.0*dwdx20;
    JB[202] = 1.0*dwdx21;
    JB[204] = -1.0*dwdx21;
    JB[205] = 1.0*dwdx21;
    JB[206] = -1.0*dwdx21;
    JB[231] = 1.0*dwdx22;
    JB[234] = -1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24;
    JB[235] = -1.0*dwdx22 - 1.0*dwdx23;
    JB[236] = -1.0*dwdx24;
    JB[238] = 1.0*dwdx24;
    JB[256] = 1.0*dwdx25;
    JB[259] = -1.0*dwdx25 + 1.0*dwdx26;
    JB[260] = -1.0*dwdx25 - 1.0*dwdx26 + 1.0*dwdx27;
    JB[261] = 1.0*dwdx27;
    JB[262] = -1.0*dwdx27;
    JB[263] = -1.0*dwdx27;
    JB[284] = 1.0*dwdx30;
    JB[285] = 1.0*dwdx28;
    JB[286] = 1.0*dwdx28 + 1.0*dwdx29 - 1.0*dwdx30;
    JB[287] = -1.0*dwdx28;
    JB[288] = -1.0*dwdx28 - 1.0*dwdx29 + 1.0*dwdx30;
    JB[293] = -1.0*dwdx29;
    JB[302] = -1.0*dwdx32;
    JB[304] = 1.0*dwdx32;
    JB[310] = 1.0*dwdx31;
    JB[311] = 1.0*dwdx31;
    JB[312] = -1.0*dwdx31 + 1.0*dwdx32;
    JB[313] = -1.0*dwdx31;
    JB[314] = -1.0*dwdx32;
    JB[334] = 1.0*dwdx35;
    JB[335] = 1.0*dwdx33;
    JB[336] = 1.0*dwdx33 + 1.0*dwdx34 - 1.0*dwdx35;
    JB[337] = -1.0*dwdx33;
    JB[338] = -1.0*dwdx33 - 1.0*dwdx34 + 1.0*dwdx35;
    JB[343] = -1.0*dwdx34;
    JB[352] = -1.0*dwdx36;
    JB[354] = 1.0*dwdx36;
    JB[362] = 1.0*dwdx36;
    JB[364] = -1.0*dwdx36 + 1.0*dwdx37;
    JB[365] = -1.0*dwdx37;
    JB[389] = 1.0*dwdx38;
    JB[390] = -1.0*dwdx38 + 1.0*dwdx39;
    JB[391] = -1.0*dwdx39;
    JB[402] = -1.0*dwdx41;
    JB[404] = 1.0*dwdx41;
    JB[415] = 1.0*dwdx40;
    JB[416] = -1.0*dwdx40 + 1.0*dwdx41;
    JB[417] = -1.0*dwdx41;
    JB[427] = -1.0*dwdx42;
    JB[429] = 1.0*dwdx42;
    JB[441] = 1.0*dwdx42;
    JB[442] = -1.0*dwdx42 + 1.0*dwdx43;
    JB[443] = -1.0*dwdx43;
    JB[461] = 1.0*dwdx44 + 3.0*dwdx45;
    JB[463] = -1.0*dwdx44 - 3.0*dwdx45;
    JB[468] = -1.0*dwdx44 + 2.0*dwdx45;
    JB[511] = 1.0*dwdx46;
    JB[513] = -1.0*dwdx46;
    JB[518] = -1.0*dwdx46;
    JB[534] = 1.0*dwdx47;
    JB[536] = -1.0*dwdx47;
    JB[538] = 1.0*dwdx47;
}