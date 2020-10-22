#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_sarma3(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 + 1.0*dwdx1;
    J[2] = -1.0*dwdx8;
    J[3] = 1.0*dwdx9;
    J[4] = 1.0*dwdx15;
    J[5] = 1.0*dwdx22;
    J[6] = 1.0*dwdx29;
    J[10] = 1.0*dwdx50;
    J[33] = 1.0*dwdx0 - 1.0*dwdx1;
    J[35] = 1.0*dwdx8;
    J[36] = -1.0*dwdx9;
    J[37] = -1.0*dwdx15;
    J[38] = -1.0*dwdx22;
    J[39] = -1.0*dwdx29;
    J[43] = -1.0*dwdx50;
    J[44] = 1.0*dwdx3;
    J[45] = 1.0*dwdx6;
    J[47] = -1.0*dwdx10 + 1.0*dwdx13;
    J[48] = -1.0*dwdx16 + 1.0*dwdx19;
    J[49] = -1.0*dwdx23 + 1.0*dwdx26;
    J[50] = 1.0*dwdx31;
    J[51] = 1.0*dwdx36;
    J[52] = 1.0*dwdx42;
    J[53] = 1.0*dwdx48;
    J[54] = 1.0*dwdx52;
    J[55] = 1.0*dwdx2 - 1.0*dwdx3;
    J[56] = 1.0*dwdx5 - 1.0*dwdx6;
    J[58] = 1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12 - 1.0*dwdx13;
    J[59] = 1.0*dwdx16 - 1.0*dwdx17 + 1.0*dwdx18 - 1.0*dwdx19;
    J[60] = 1.0*dwdx23 - 1.0*dwdx24 + 1.0*dwdx25 - 1.0*dwdx26;
    J[61] = 1.0*dwdx30 - 1.0*dwdx31;
    J[62] = 1.0*dwdx35 - 1.0*dwdx36;
    J[63] = 1.0*dwdx41 - 1.0*dwdx42;
    J[64] = 1.0*dwdx47 - 1.0*dwdx48;
    J[65] = 1.0*dwdx51 - 1.0*dwdx52;
    J[66] = -1.0*dwdx2;
    J[67] = -1.0*dwdx5;
    J[69] = 1.0*dwdx11 - 1.0*dwdx12;
    J[70] = 1.0*dwdx17 - 1.0*dwdx18;
    J[71] = 1.0*dwdx24 - 1.0*dwdx25;
    J[72] = -1.0*dwdx30;
    J[73] = -1.0*dwdx35;
    J[74] = -1.0*dwdx41;
    J[75] = -1.0*dwdx47;
    J[76] = -1.0*dwdx51;
    J[78] = 1.0*dwdx4;
    J[81] = 1.0*dwdx14;
    J[82] = 1.0*dwdx21;
    J[83] = 1.0*dwdx28 - 1.0*dwdx32;
    J[84] = 1.0*dwdx34 - 1.0*dwdx37;
    J[85] = 1.0*dwdx40 - 1.0*dwdx43;
    J[86] = 1.0*dwdx46;
    J[89] = -1.0*dwdx4 + 1.0*dwdx7;
    J[92] = -1.0*dwdx14 + 1.0*dwdx20;
    J[93] = -1.0*dwdx21 + 1.0*dwdx27;
    J[94] = -1.0*dwdx28 + 1.0*dwdx32 - 1.0*dwdx33;
    J[95] = -1.0*dwdx34 + 1.0*dwdx37 - 1.0*dwdx38 + 1.0*dwdx39;
    J[96] = -1.0*dwdx40 + 1.0*dwdx43 - 1.0*dwdx44 + 1.0*dwdx45;
    J[97] = -1.0*dwdx46 + 1.0*dwdx49;
    J[98] = 1.0*dwdx53;
    J[100] = -1.0*dwdx7;
    J[103] = -1.0*dwdx20;
    J[104] = -1.0*dwdx27;
    J[105] = 1.0*dwdx33;
    J[106] = 1.0*dwdx38 - 1.0*dwdx39;
    J[107] = 1.0*dwdx44 - 1.0*dwdx45;
    J[108] = -1.0*dwdx49;
    J[109] = -1.0*dwdx53;
}