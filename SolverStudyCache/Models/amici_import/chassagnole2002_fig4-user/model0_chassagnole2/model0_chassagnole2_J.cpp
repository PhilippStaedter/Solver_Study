#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_chassagnole2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    J[3] = 1.0*dwdx13;
    J[6] = 1.0*dwdx26 - 1.0*dwdx30;
    J[19] = -1.0*dwdx4 - 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7;
    J[20] = 1.0*dwdx10 - 1.0*dwdx11;
    J[24] = 1.0*dwdx29 - 1.0*dwdx32;
    J[26] = -1.0*dwdx35;
    J[34] = 1.0*dwdx68;
    J[35] = -1.0*dwdx72;
    J[37] = 1.0*dwdx6 + 1.0*dwdx7;
    J[38] = 1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx8 + 1.0*dwdx9;
    J[41] = 1.0*dwdx23;
    J[42] = 1.0*dwdx29 + 1.0*dwdx32;
    J[44] = -1.0*dwdx38;
    J[45] = 1.0*dwdx45;
    J[52] = 1.0*dwdx68;
    J[53] = 1.0*dwdx72;
    J[54] = -1.0*dwdx0;
    J[56] = 1.0*dwdx8;
    J[57] = -1.0*dwdx13 - 1.0*dwdx16;
    J[60] = -1.0*dwdx26;
    J[62] = 1.0*dwdx38;
    J[75] = -1.0*dwdx14;
    J[76] = -1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
    J[77] = 1.0*dwdx24;
    J[92] = -1.0*dwdx9;
    J[94] = -1.0*dwdx20;
    J[95] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 + 65.0*dwdx25;
    J[97] = 65.0*dwdx34;
    J[98] = 65.0*dwdx40;
    J[99] = -1.0*dwdx45;
    J[103] = 65.0*dwdx57;
    J[108] = 1.0*dwdx0 + 1.0*dwdx3;
    J[109] = -1.0*dwdx6 + 1.0*dwdx7;
    J[110] = -1.0*dwdx10 + 1.0*dwdx11;
    J[111] = 1.0*dwdx13;
    J[114] = 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31 + 1.0*dwdx32;
    J[120] = -1.0*dwdx53;
    J[122] = 1.0*dwdx63;
    J[124] = -1.0*dwdx68 + 1.0*dwdx69;
    J[125] = 1.0*dwdx71 + 1.0*dwdx72;
    J[131] = -1.0*dwdx25;
    J[133] = 1.0*dwdx33 - 1.0*dwdx34;
    J[134] = -1.0*dwdx40;
    J[139] = -1.0*dwdx57;
    J[145] = -1.0*dwdx4;
    J[147] = -1.0*dwdx15 - 1.0*dwdx17;
    J[149] = -65.0*dwdx25;
    J[151] = -65.0*dwdx34;
    J[152] = -1.0*dwdx35 + 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx39 - 65.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42;
    J[154] = 1.0*dwdx46;
    J[157] = -65.0*dwdx57;
    J[167] = 1.0*dwdx22;
    J[171] = -1.0*dwdx43 - 1.0*dwdx44;
    J[188] = -1.0*dwdx36;
    J[190] = -1.0*dwdx46 - 1.0*dwdx47 + 1.0*dwdx48;
    J[191] = 1.0*dwdx51;
    J[208] = -1.0*dwdx48;
    J[209] = -1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    J[210] = 1.0*dwdx54;
    J[222] = 1.0*dwdx28;
    J[227] = -1.0*dwdx50;
    J[228] = 1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55;
    J[237] = 1.0*dwdx15;
    J[239] = 65.0*dwdx25;
    J[241] = 65.0*dwdx34;
    J[242] = 1.0*dwdx39 + 65.0*dwdx40;
    J[247] = -1.0*dwdx56 + 65.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    J[258] = -1.0*dwdx31;
    J[266] = -1.0*dwdx60 + 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63;
    J[267] = 1.0*dwdx64;
    J[268] = -1.0*dwdx69;
    J[269] = -1.0*dwdx71;
    J[279] = 1.0*dwdx44;
    J[284] = -1.0*dwdx61;
    J[285] = -1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66;
    J[287] = -1.0*dwdx70;
    J[289] = -1.0*dwdx6;
    J[290] = -1.0*dwdx10;
    J[294] = -1.0*dwdx29 + 1.0*dwdx31;
    J[302] = 1.0*dwdx63;
    J[304] = -1.0*dwdx67 - 1.0*dwdx68 + 1.0*dwdx69;
    J[305] = 1.0*dwdx71;
    J[307] = -1.0*dwdx7;
    J[308] = -1.0*dwdx11;
    J[312] = -1.0*dwdx31 - 1.0*dwdx32;
    J[320] = -1.0*dwdx63;
    J[321] = 1.0*dwdx66;
    J[322] = -1.0*dwdx69;
    J[323] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73;
}