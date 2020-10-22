#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Leber2015(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx3 - 1.0*dwdx7;
    J[1] = -1.0*dwdx10;
    J[2] = -1.0*dwdx13;
    J[5] = -1.0*dwdx19 + 1.0*dwdx20;
    J[6] = -1.0*dwdx24;
    J[7] = -1.0*dwdx29;
    J[10] = -1.0*dwdx35;
    J[11] = -1.0*dwdx39;
    J[24] = -1.0*dwdx9;
    J[25] = -1.0*dwdx12;
    J[29] = -1.0*dwdx26;
    J[33] = -1.0*dwdx34;
    J[47] = 1.0*dwdx9;
    J[48] = 1.0*dwdx12 - 1.0*dwdx15;
    J[52] = 1.0*dwdx26;
    J[56] = 1.0*dwdx34;
    J[69] = 1.0*dwdx7;
    J[70] = 1.0*dwdx10;
    J[71] = 1.0*dwdx13;
    J[72] = -1.0*dwdx16;
    J[76] = 1.0*dwdx29;
    J[79] = 1.0*dwdx35;
    J[95] = 1.0*dwdx16;
    J[96] = -1.0*dwdx17 - 1.0*dwdx18;
    J[120] = -1.0*dwdx21;
    J[125] = -1.0*dwdx36;
    J[129] = -1.0*dwdx45;
    J[138] = 1.0*dwdx2;
    J[139] = -1.0*dwdx8;
    J[144] = -1.0*dwdx22;
    J[146] = 1.0*dwdx31;
    J[152] = 1.0*dwdx44;
    J[153] = 1.0*dwdx48;
    J[155] = 1.0*dwdx54;
    J[161] = -0.25*dwdx5;
    J[167] = -0.25*dwdx23;
    J[168] = -0.25*dwdx27 - 0.25*dwdx28;
    J[169] = 0.25*dwdx32;
    J[172] = -0.25*dwdx38;
    J[176] = -0.25*dwdx47;
    J[190] = 0.25*dwdx23 + 0.25*dwdx25;
    J[191] = 0.25*dwdx27;
    J[192] = -0.25*dwdx32;
    J[194] = 0.25*dwdx33 + 0.25*dwdx37;
    J[195] = 0.25*dwdx38 + 0.25*dwdx40;
    J[199] = 0.25*dwdx47 + 0.25*dwdx50;
    J[230] = 0.25*dwdx5;
    J[236] = -0.25*dwdx25;
    J[237] = 0.25*dwdx28;
    J[240] = -0.25*dwdx33 - 0.25*dwdx37;
    J[241] = -0.25*dwdx40;
    J[245] = -0.25*dwdx50;
    J[253] = 0.25*dwdx6;
    J[264] = -0.25*dwdx41;
    J[266] = 0.25*dwdx43;
    J[268] = 0.25*dwdx51;
    J[270] = 0.25*dwdx56;
    J[276] = 14.285714285714285*dwdx0;
    J[288] = -14.285714285714285*dwdx42;
    J[345] = -14.285714285714285*dwdx4;
    J[360] = -14.285714285714285*dwdx46 - 14.285714285714285*dwdx49;
    J[362] = -14.285714285714285*dwdx55;
    J[366] = 14.285714285714285*dwdx61;
    J[384] = -14.285714285714285*dwdx52;
    J[390] = 14.285714285714285*dwdx62;
    J[391] = 14.285714285714285*dwdx4;
    J[406] = 14.285714285714285*dwdx49;
    J[408] = -14.285714285714285*dwdx53 + 14.285714285714285*dwdx55;
    J[410] = 14.285714285714285*dwdx60;
    J[415] = -1.0*dwdx11;
    J[416] = -1.0*dwdx14;
    J[421] = -1.0*dwdx30;
    J[426] = 1.0*dwdx42;
    J[432] = -1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    J[441] = 1.0*dwdx18;
    J[456] = -1.0*dwdx60;
    J[501] = 1.0*dwdx58;
    J[504] = -1.0*dwdx61;
    J[507] = 1.0*dwdx11;
    J[508] = 1.0*dwdx14;
    J[513] = 1.0*dwdx30;
    J[524] = 1.0*dwdx59;
    J[528] = -1.0*dwdx62;
}