#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_chassagnole2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[3] = 1.0*dwdx0;
    JB[6] = -1.0*dwdx0 - 1.0*dwdx3;
    JB[19] = 1.0*dwdx4 + 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7;
    JB[20] = -1.0*dwdx6 - 1.0*dwdx7;
    JB[24] = 1.0*dwdx6 - 1.0*dwdx7;
    JB[26] = 1.0*dwdx4;
    JB[34] = 1.0*dwdx6;
    JB[35] = 1.0*dwdx7;
    JB[37] = -1.0*dwdx10 + 1.0*dwdx11;
    JB[38] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx8 - 1.0*dwdx9;
    JB[39] = -1.0*dwdx8;
    JB[41] = 1.0*dwdx9;
    JB[42] = 1.0*dwdx10 - 1.0*dwdx11;
    JB[52] = 1.0*dwdx10;
    JB[53] = 1.0*dwdx11;
    JB[54] = -1.0*dwdx13;
    JB[57] = 1.0*dwdx13 + 1.0*dwdx16;
    JB[58] = 1.0*dwdx14;
    JB[60] = -1.0*dwdx13;
    JB[62] = 1.0*dwdx15 + 1.0*dwdx17;
    JB[67] = -1.0*dwdx15;
    JB[76] = 1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20;
    JB[77] = 1.0*dwdx20;
    JB[92] = -1.0*dwdx23;
    JB[94] = -1.0*dwdx24;
    JB[95] = 1.0*dwdx21 + 1.0*dwdx22 + 1.0*dwdx23 + 1.0*dwdx24 - 65.0*dwdx25;
    JB[97] = 1.0*dwdx25;
    JB[98] = 65.0*dwdx25;
    JB[99] = -1.0*dwdx22;
    JB[103] = -65.0*dwdx25;
    JB[108] = -1.0*dwdx26 + 1.0*dwdx30;
    JB[109] = -1.0*dwdx29 + 1.0*dwdx32;
    JB[110] = -1.0*dwdx29 - 1.0*dwdx32;
    JB[111] = 1.0*dwdx26;
    JB[114] = -1.0*dwdx26 + 1.0*dwdx27 + 1.0*dwdx28 + 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32;
    JB[120] = -1.0*dwdx28;
    JB[122] = 1.0*dwdx31;
    JB[124] = 1.0*dwdx29 - 1.0*dwdx31;
    JB[125] = 1.0*dwdx31 + 1.0*dwdx32;
    JB[131] = -65.0*dwdx34;
    JB[133] = -1.0*dwdx33 + 1.0*dwdx34;
    JB[134] = 65.0*dwdx34;
    JB[139] = -65.0*dwdx34;
    JB[145] = 1.0*dwdx35;
    JB[146] = 1.0*dwdx38;
    JB[147] = -1.0*dwdx38;
    JB[149] = -65.0*dwdx40;
    JB[151] = 1.0*dwdx40;
    JB[152] = 1.0*dwdx35 - 1.0*dwdx36 + 1.0*dwdx37 + 1.0*dwdx39 + 65.0*dwdx40 + 1.0*dwdx41 + 1.0*dwdx42;
    JB[154] = 1.0*dwdx36;
    JB[157] = -1.0*dwdx39 - 65.0*dwdx40;
    JB[164] = -1.0*dwdx45;
    JB[167] = 1.0*dwdx45;
    JB[171] = 1.0*dwdx43 + 1.0*dwdx44;
    JB[177] = -1.0*dwdx44;
    JB[188] = -1.0*dwdx46;
    JB[190] = 1.0*dwdx46 + 1.0*dwdx47 - 1.0*dwdx48;
    JB[191] = 1.0*dwdx48;
    JB[208] = -1.0*dwdx51;
    JB[209] = 1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51 + 1.0*dwdx52;
    JB[210] = 1.0*dwdx50;
    JB[222] = 1.0*dwdx53;
    JB[227] = -1.0*dwdx54;
    JB[228] = -1.0*dwdx53 + 1.0*dwdx54 + 1.0*dwdx55;
    JB[239] = -65.0*dwdx57;
    JB[241] = 1.0*dwdx57;
    JB[242] = 65.0*dwdx57;
    JB[247] = 1.0*dwdx56 - 65.0*dwdx57 + 1.0*dwdx58 + 1.0*dwdx59;
    JB[258] = -1.0*dwdx63;
    JB[266] = 1.0*dwdx60 - 1.0*dwdx61 + 1.0*dwdx62 + 1.0*dwdx63;
    JB[267] = 1.0*dwdx61;
    JB[268] = -1.0*dwdx63;
    JB[269] = 1.0*dwdx63;
    JB[284] = -1.0*dwdx64;
    JB[285] = 1.0*dwdx64 + 1.0*dwdx65 + 1.0*dwdx66;
    JB[287] = -1.0*dwdx66;
    JB[289] = -1.0*dwdx68;
    JB[290] = -1.0*dwdx68;
    JB[294] = 1.0*dwdx68 - 1.0*dwdx69;
    JB[302] = 1.0*dwdx69;
    JB[304] = 1.0*dwdx67 + 1.0*dwdx68 - 1.0*dwdx69;
    JB[305] = 1.0*dwdx69;
    JB[307] = 1.0*dwdx72;
    JB[308] = -1.0*dwdx72;
    JB[312] = -1.0*dwdx71 - 1.0*dwdx72;
    JB[320] = 1.0*dwdx71;
    JB[321] = 1.0*dwdx70;
    JB[322] = -1.0*dwdx71;
    JB[323] = -1.0*dwdx70 + 1.0*dwdx71 + 1.0*dwdx72 + 1.0*dwdx73;
}