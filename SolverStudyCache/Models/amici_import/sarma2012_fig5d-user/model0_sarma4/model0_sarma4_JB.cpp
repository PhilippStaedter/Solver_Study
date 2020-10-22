#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_sarma4(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 - 1.0*dwdx1;
    JB[3] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[4] = -1.0*dwdx3;
    JB[5] = -1.0*dwdx2 + 1.0*dwdx3;
    JB[6] = 1.0*dwdx2;
    JB[15] = -1.0*dwdx6;
    JB[16] = -1.0*dwdx5 + 1.0*dwdx6;
    JB[17] = 1.0*dwdx5;
    JB[18] = -1.0*dwdx4;
    JB[19] = 1.0*dwdx4 - 1.0*dwdx7;
    JB[20] = 1.0*dwdx7;
    JB[22] = 1.0*dwdx8;
    JB[25] = -1.0*dwdx8;
    JB[33] = -1.0*dwdx9;
    JB[36] = 1.0*dwdx9;
    JB[37] = 1.0*dwdx10 - 1.0*dwdx13;
    JB[38] = -1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx12 + 1.0*dwdx13;
    JB[39] = -1.0*dwdx11 + 1.0*dwdx12;
    JB[44] = -1.0*dwdx15;
    JB[47] = 1.0*dwdx15;
    JB[48] = 1.0*dwdx16 - 1.0*dwdx19;
    JB[49] = -1.0*dwdx16 + 1.0*dwdx17 - 1.0*dwdx18 + 1.0*dwdx19;
    JB[50] = -1.0*dwdx17 + 1.0*dwdx18;
    JB[51] = -1.0*dwdx14;
    JB[52] = 1.0*dwdx14 - 1.0*dwdx20;
    JB[53] = 1.0*dwdx20;
    JB[55] = -1.0*dwdx22;
    JB[58] = 1.0*dwdx22;
    JB[59] = 1.0*dwdx23 - 1.0*dwdx26;
    JB[60] = -1.0*dwdx23 + 1.0*dwdx24 - 1.0*dwdx25 + 1.0*dwdx26;
    JB[61] = -1.0*dwdx24 + 1.0*dwdx25;
    JB[62] = -1.0*dwdx21;
    JB[63] = 1.0*dwdx21 - 1.0*dwdx27;
    JB[64] = 1.0*dwdx27;
    JB[66] = -1.0*dwdx29;
    JB[69] = 1.0*dwdx29;
    JB[70] = -1.0*dwdx31;
    JB[71] = -1.0*dwdx30 + 1.0*dwdx31;
    JB[72] = 1.0*dwdx30;
    JB[73] = -1.0*dwdx28 + 1.0*dwdx32;
    JB[74] = 1.0*dwdx28 - 1.0*dwdx32 + 1.0*dwdx33;
    JB[75] = -1.0*dwdx33;
    JB[81] = -1.0*dwdx36;
    JB[82] = -1.0*dwdx35 + 1.0*dwdx36;
    JB[83] = 1.0*dwdx35;
    JB[84] = -1.0*dwdx34 + 1.0*dwdx37;
    JB[85] = 1.0*dwdx34 - 1.0*dwdx37 + 1.0*dwdx38 - 1.0*dwdx39;
    JB[86] = -1.0*dwdx38 + 1.0*dwdx39;
    JB[92] = -1.0*dwdx42;
    JB[93] = -1.0*dwdx41 + 1.0*dwdx42;
    JB[94] = 1.0*dwdx41;
    JB[95] = -1.0*dwdx40 + 1.0*dwdx43;
    JB[96] = 1.0*dwdx40 - 1.0*dwdx43 + 1.0*dwdx44 - 1.0*dwdx45;
    JB[97] = -1.0*dwdx44 + 1.0*dwdx45;
    JB[103] = -1.0*dwdx48;
    JB[104] = -1.0*dwdx47 + 1.0*dwdx48;
    JB[105] = 1.0*dwdx47;
    JB[106] = -1.0*dwdx46;
    JB[107] = 1.0*dwdx46 - 1.0*dwdx49;
    JB[108] = 1.0*dwdx49;
    JB[110] = -1.0*dwdx50;
    JB[113] = 1.0*dwdx50;
    JB[114] = -1.0*dwdx52;
    JB[115] = -1.0*dwdx51 + 1.0*dwdx52;
    JB[116] = 1.0*dwdx51;
    JB[118] = -1.0*dwdx53;
    JB[119] = 1.0*dwdx53;
}