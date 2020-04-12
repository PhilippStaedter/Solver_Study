#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_Leber2015(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 - 1.0*dwdx3 + 1.0*dwdx7;
    JB[3] = -1.0*dwdx7;
    JB[6] = -1.0*dwdx2;
    JB[7] = 0.25*dwdx5;
    JB[10] = -0.25*dwdx5;
    JB[11] = -0.25*dwdx6;
    JB[12] = -14.285714285714285*dwdx0;
    JB[15] = 14.285714285714285*dwdx4;
    JB[17] = -14.285714285714285*dwdx4;
    JB[23] = 1.0*dwdx10;
    JB[24] = 1.0*dwdx9;
    JB[25] = -1.0*dwdx9;
    JB[26] = -1.0*dwdx10;
    JB[29] = 1.0*dwdx8;
    JB[41] = 1.0*dwdx11;
    JB[45] = -1.0*dwdx11;
    JB[46] = 1.0*dwdx13;
    JB[47] = 1.0*dwdx12;
    JB[48] = -1.0*dwdx12 + 1.0*dwdx15;
    JB[49] = -1.0*dwdx13;
    JB[64] = 1.0*dwdx14;
    JB[68] = -1.0*dwdx14;
    JB[72] = 1.0*dwdx16;
    JB[73] = -1.0*dwdx16;
    JB[96] = 1.0*dwdx17 + 1.0*dwdx18;
    JB[111] = -1.0*dwdx18;
    JB[115] = 1.0*dwdx19 - 1.0*dwdx20;
    JB[120] = 1.0*dwdx21;
    JB[138] = 1.0*dwdx24;
    JB[139] = 1.0*dwdx26;
    JB[140] = -1.0*dwdx26;
    JB[144] = 1.0*dwdx22;
    JB[145] = 0.25*dwdx23;
    JB[146] = -0.25*dwdx23 - 0.25*dwdx25;
    JB[148] = 0.25*dwdx25;
    JB[161] = 1.0*dwdx29;
    JB[164] = -1.0*dwdx29;
    JB[168] = 0.25*dwdx27 + 0.25*dwdx28;
    JB[169] = -0.25*dwdx27;
    JB[171] = -0.25*dwdx28;
    JB[179] = 1.0*dwdx30;
    JB[183] = -1.0*dwdx30;
    JB[190] = -1.0*dwdx31;
    JB[191] = -0.25*dwdx32;
    JB[192] = 0.25*dwdx32;
    JB[230] = 1.0*dwdx35;
    JB[231] = 1.0*dwdx34;
    JB[232] = -1.0*dwdx34;
    JB[233] = -1.0*dwdx35;
    JB[235] = 1.0*dwdx36;
    JB[238] = -0.25*dwdx33 - 0.25*dwdx37;
    JB[240] = 0.25*dwdx33 + 0.25*dwdx37;
    JB[253] = 1.0*dwdx39;
    JB[260] = 0.25*dwdx38;
    JB[261] = -0.25*dwdx38 - 0.25*dwdx40;
    JB[263] = 0.25*dwdx40;
    JB[264] = 0.25*dwdx41;
    JB[288] = 14.285714285714285*dwdx42;
    JB[294] = -1.0*dwdx42;
    JB[310] = -0.25*dwdx43;
    JB[327] = 1.0*dwdx45;
    JB[328] = -1.0*dwdx44;
    JB[351] = -1.0*dwdx48;
    JB[352] = 0.25*dwdx47;
    JB[353] = -0.25*dwdx47 - 0.25*dwdx50;
    JB[355] = 0.25*dwdx50;
    JB[356] = -0.25*dwdx51;
    JB[360] = 14.285714285714285*dwdx46 + 14.285714285714285*dwdx49;
    JB[362] = -14.285714285714285*dwdx49;
    JB[384] = 14.285714285714285*dwdx52;
    JB[397] = -1.0*dwdx54;
    JB[402] = -0.25*dwdx56;
    JB[406] = 14.285714285714285*dwdx55;
    JB[408] = 14.285714285714285*dwdx53 - 14.285714285714285*dwdx55;
    JB[432] = 1.0*dwdx57 + 1.0*dwdx58 + 1.0*dwdx59;
    JB[435] = -1.0*dwdx58;
    JB[436] = -1.0*dwdx59;
    JB[454] = -14.285714285714285*dwdx60;
    JB[456] = 1.0*dwdx60;
    JB[498] = -14.285714285714285*dwdx61;
    JB[504] = 1.0*dwdx61;
    JB[522] = -14.285714285714285*dwdx62;
    JB[528] = 1.0*dwdx62;
}