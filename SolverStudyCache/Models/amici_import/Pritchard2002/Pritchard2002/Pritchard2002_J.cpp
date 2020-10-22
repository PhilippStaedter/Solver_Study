#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Pritchard2002(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[25] = 1.0*dwdx0;
    J[26] = 1.0*dwdx1 - 1.0*dwdx2;
    J[27] = -1.0*dwdx3;
    J[28] = -1.0*dwdx9;
    J[29] = -1.0*dwdx11;
    J[51] = -1.0*dwdx2;
    J[52] = -1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7 + 1.0*dwdx8;
    J[53] = -1.0*dwdx9;
    J[54] = -1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14;
    J[55] = -1.0*dwdx16;
    J[56] = -1.0*dwdx17;
    J[57] = -1.0*dwdx19 + 1.0*dwdx20;
    J[58] = -1.0*dwdx21;
    J[62] = 1.0*dwdx32;
    J[64] = 1.0*dwdx36;
    J[66] = 1.0*dwdx41;
    J[67] = 1.0*dwdx42;
    J[76] = 1.0*dwdx2;
    J[77] = 1.0*dwdx3;
    J[78] = -1.0*dwdx10 + 1.0*dwdx9;
    J[79] = 1.0*dwdx11;
    J[80] = -1.0*dwdx15;
    J[101] = 1.0*dwdx2;
    J[102] = 1.0*dwdx3 + 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7 - 2.0*dwdx8;
    J[103] = 1.0*dwdx9;
    J[104] = 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 2.0*dwdx14;
    J[105] = 1.0*dwdx16;
    J[106] = 1.0*dwdx17;
    J[107] = 1.0*dwdx19 - 2.0*dwdx20;
    J[108] = 1.0*dwdx21;
    J[112] = -1.0*dwdx32;
    J[114] = -1.0*dwdx36;
    J[116] = -1.0*dwdx41;
    J[117] = -1.0*dwdx42;
    J[127] = -1.0*dwdx4;
    J[128] = 1.0*dwdx10;
    J[130] = 1.0*dwdx15 - 1.0*dwdx16;
    J[131] = -1.0*dwdx17;
    J[132] = -1.0*dwdx19;
    J[133] = -1.0*dwdx21;
    J[152] = 1.0*dwdx4;
    J[155] = 1.0*dwdx16;
    J[156] = 1.0*dwdx17 - 1.0*dwdx18;
    J[157] = 1.0*dwdx19;
    J[158] = 1.0*dwdx21;
    J[159] = -1.0*dwdx22;
    J[160] = -1.0*dwdx25;
    J[177] = 1.0*dwdx8;
    J[179] = 1.0*dwdx14;
    J[182] = 1.0*dwdx20;
    J[231] = 1.0*dwdx18;
    J[234] = 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    J[235] = 1.0*dwdx25 - 1.0*dwdx26;
    J[236] = -1.0*dwdx30;
    J[238] = -1.0*dwdx35;
    J[246] = -1.0*dwdx47;
    J[256] = 1.0*dwdx18;
    J[259] = 1.0*dwdx22 + 1.0*dwdx23;
    J[260] = 1.0*dwdx25 + 1.0*dwdx26 - 1.0*dwdx27;
    J[261] = -1.0*dwdx28;
    J[262] = -1.0*dwdx31;
    J[263] = -1.0*dwdx33;
    J[284] = 1.0*dwdx24;
    J[285] = -1.0*dwdx27;
    J[286] = -1.0*dwdx28 - 1.0*dwdx29 + 1.0*dwdx30;
    J[287] = -1.0*dwdx31;
    J[288] = -1.0*dwdx33 - 1.0*dwdx34 + 1.0*dwdx35;
    J[293] = -1.0*dwdx44 - 3.0*dwdx45;
    J[295] = -1.0*dwdx46;
    J[296] = 1.0*dwdx47;
    J[302] = -1.0*dwdx5;
    J[304] = -1.0*dwdx12;
    J[310] = 1.0*dwdx27;
    J[311] = 1.0*dwdx28;
    J[312] = 1.0*dwdx31 - 1.0*dwdx32;
    J[313] = 1.0*dwdx33;
    J[314] = -1.0*dwdx36;
    J[334] = -1.0*dwdx24;
    J[335] = 1.0*dwdx27;
    J[336] = 1.0*dwdx28 + 1.0*dwdx29 - 1.0*dwdx30;
    J[337] = 1.0*dwdx31;
    J[338] = 1.0*dwdx33 + 1.0*dwdx34 - 1.0*dwdx35;
    J[343] = 1.0*dwdx44 + 3.0*dwdx45;
    J[345] = 1.0*dwdx46;
    J[346] = -1.0*dwdx47;
    J[352] = 1.0*dwdx5;
    J[354] = 1.0*dwdx12;
    J[362] = 1.0*dwdx32;
    J[364] = 1.0*dwdx36 - 1.0*dwdx37;
    J[365] = -1.0*dwdx38;
    J[389] = 1.0*dwdx37;
    J[390] = 1.0*dwdx38 - 1.0*dwdx39;
    J[391] = -1.0*dwdx40;
    J[402] = -1.0*dwdx6;
    J[404] = -1.0*dwdx13;
    J[415] = 1.0*dwdx39;
    J[416] = 1.0*dwdx40 - 1.0*dwdx41;
    J[417] = -1.0*dwdx42;
    J[427] = 1.0*dwdx6;
    J[429] = 1.0*dwdx13;
    J[441] = 1.0*dwdx41;
    J[442] = 1.0*dwdx42 - 1.0*dwdx43;
    J[461] = 1.0*dwdx29;
    J[463] = 1.0*dwdx34;
    J[467] = 1.0*dwdx43;
    J[468] = 1.0*dwdx44 - 2.0*dwdx45;
    J[470] = 1.0*dwdx46;
}