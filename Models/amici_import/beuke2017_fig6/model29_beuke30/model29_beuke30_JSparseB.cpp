#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model29_beuke30(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 76923.076923076922*dwdx0 + 76923.076923076922*dwdx1 + 76923.076923076922*dwdx2 + 76923.076923076922*dwdx3;
    JSparseB[1] = -76923.076923076922*dwdx5;
    JSparseB[2] = 76923.076923076922*dwdx47;
    JSparseB[3] = 76923.076923076922*dwdx55;
    JSparseB[4] = 76923.076923076922*dwdx59;
    JSparseB[5] = 76923.076923076922*dwdx64 - 76923.076923076922*dwdx66;
    JSparseB[6] = 76923.076923076922*dwdx4;
    JSparseB[7] = -76923.076923076922*dwdx7;
    JSparseB[8] = 1250000.0*dwdx6 + 1250000.0*dwdx7;
    JSparseB[9] = -1250000.0*dwdx9;
    JSparseB[10] = 1250000.0*dwdx8;
    JSparseB[11] = -1250000.0*dwdx42;
    JSparseB[12] = -1250000.0*dwdx52;
    JSparseB[13] = -1250000.0*dwdx76;
    JSparseB[14] = -1250000.0*dwdx8;
    JSparseB[15] = 1250000.0*dwdx9;
    JSparseB[16] = 76923.076923076922*dwdx10 + 76923.076923076922*dwdx11;
    JSparseB[17] = -76923.076923076922*dwdx12;
    JSparseB[18] = 76923.076923076922*dwdx25;
    JSparseB[19] = -76923.076923076922*dwdx10 - 76923.076923076922*dwdx11;
    JSparseB[20] = 76923.076923076922*dwdx12;
    JSparseB[21] = -76923.076923076922*dwdx25;
    JSparseB[22] = 386100.38610038609*dwdx13 + 386100.38610038609*dwdx14;
    JSparseB[23] = 1250000.0*dwdx17;
    JSparseB[24] = -1250000.0*dwdx18;
    JSparseB[25] = 1250000.0*dwdx38;
    JSparseB[26] = -1250000.0*dwdx17;
    JSparseB[27] = 1250000.0*dwdx18;
    JSparseB[28] = -1250000.0*dwdx38;
    JSparseB[29] = -386100.38610038609*dwdx21 + 386100.38610038609*dwdx22;
    JSparseB[30] = -386100.38610038609*dwdx23;
    JSparseB[31] = 386100.38610038609*dwdx26;
    JSparseB[32] = 386100.38610038609*dwdx62;
    JSparseB[33] = 76923.076923076922*dwdx21;
    JSparseB[34] = 76923.076923076922*dwdx23;
    JSparseB[35] = -76923.076923076922*dwdx29;
    JSparseB[36] = -386100.38610038609*dwdx22;
    JSparseB[37] = -386100.38610038609*dwdx26 + 386100.38610038609*dwdx27;
    JSparseB[38] = -386100.38610038609*dwdx62;
    JSparseB[39] = -76923.076923076922*dwdx27;
    JSparseB[40] = 76923.076923076922*dwdx29;
    JSparseB[41] = 2127659.5744680851*dwdx31;
    JSparseB[42] = -2127659.5744680851*dwdx34;
    JSparseB[43] = 2857142.8571428573*dwdx33;
    JSparseB[44] = -2857142.8571428573*dwdx35;
    JSparseB[45] = -4255319.1489361702*dwdx15;
    JSparseB[46] = 2127659.5744680851*dwdx34;
    JSparseB[47] = -5714285.7142857146*dwdx16;
    JSparseB[48] = 2857142.8571428573*dwdx35;
    JSparseB[49] = 76923.076923076922*dwdx28;
    JSparseB[50] = 76923.076923076922*dwdx36 + 76923.076923076922*dwdx37;
    JSparseB[51] = -76923.076923076922*dwdx39;
    JSparseB[52] = -76923.076923076922*dwdx28;
    JSparseB[53] = -76923.076923076922*dwdx36 - 76923.076923076922*dwdx37;
    JSparseB[54] = 76923.076923076922*dwdx39;
    JSparseB[55] = -1250000.0*dwdx19 - 1250000.0*dwdx20;
    JSparseB[56] = 1250000.0*dwdx40 + 1250000.0*dwdx41;
    JSparseB[57] = -1250000.0*dwdx51;
    JSparseB[58] = -1250000.0*dwdx75;
    JSparseB[59] = 1250000.0*dwdx43 + 1250000.0*dwdx44 + 1250000.0*dwdx45 - 1250000.0*dwdx46;
    JSparseB[60] = -1250000.0*dwdx53;
    JSparseB[61] = -1250000.0*dwdx56;
    JSparseB[62] = 1250000.0*dwdx67;
    JSparseB[63] = 76923.076923076922*dwdx0;
    JSparseB[64] = 76923.076923076922*dwdx47 + 76923.076923076922*dwdx48;
    JSparseB[65] = -76923.076923076922*dwdx49;
    JSparseB[66] = 76923.076923076922*dwdx50;
    JSparseB[67] = 76923.076923076922*dwdx64 - 76923.076923076922*dwdx65;
    JSparseB[68] = -76923.076923076922*dwdx70;
    JSparseB[69] = 1250000.0*dwdx20;
    JSparseB[70] = -1250000.0*dwdx41;
    JSparseB[71] = -1250000.0*dwdx43 + 1250000.0*dwdx46;
    JSparseB[72] = -1250000.0*dwdx48;
    JSparseB[73] = -1250000.0*dwdx50 + 1250000.0*dwdx51 + 1250000.0*dwdx53;
    JSparseB[74] = 1250000.0*dwdx56;
    JSparseB[75] = -1250000.0*dwdx74;
    JSparseB[76] = -1250000.0*dwdx2;
    JSparseB[77] = -1250000.0*dwdx44 + 1250000.0*dwdx46;
    JSparseB[78] = 1250000.0*dwdx53;
    JSparseB[79] = 1250000.0*dwdx54 - 1250000.0*dwdx55 + 1250000.0*dwdx56;
    JSparseB[80] = 76923.076923076922*dwdx24;
    JSparseB[81] = 76923.076923076922*dwdx57 + 76923.076923076922*dwdx58;
    JSparseB[82] = -76923.076923076922*dwdx61;
    JSparseB[83] = -76923.076923076922*dwdx24;
    JSparseB[84] = -76923.076923076922*dwdx57 - 76923.076923076922*dwdx58;
    JSparseB[85] = 76923.076923076922*dwdx61;
    JSparseB[86] = 386100.38610038609*dwdx22;
    JSparseB[87] = 386100.38610038609*dwdx26;
    JSparseB[88] = -386100.38610038609*dwdx30;
    JSparseB[89] = -386100.38610038609*dwdx32;
    JSparseB[90] = 386100.38610038609*dwdx62 + 386100.38610038609*dwdx63;
    JSparseB[91] = -76923.076923076922*dwdx0;
    JSparseB[92] = -76923.076923076922*dwdx45;
    JSparseB[93] = -76923.076923076922*dwdx47;
    JSparseB[94] = 76923.076923076922*dwdx60;
    JSparseB[95] = -76923.076923076922*dwdx64 + 76923.076923076922*dwdx65 + 76923.076923076922*dwdx66 - 76923.076923076922*dwdx67 + 76923.076923076922*dwdx68;
    JSparseB[96] = -76923.076923076922*dwdx60;
    JSparseB[97] = -76923.076923076922*dwdx68;
    JSparseB[98] = 76923.076923076922*dwdx69 + 76923.076923076922*dwdx70 + 76923.076923076922*dwdx71;
    JSparseB[99] = 76923.076923076922*dwdx77;
    JSparseB[100] = -76923.076923076922*dwdx3;
    JSparseB[101] = -76923.076923076922*dwdx59 - 76923.076923076922*dwdx60;
    JSparseB[102] = -76923.076923076922*dwdx68;
    JSparseB[103] = 76923.076923076922*dwdx72;
    JSparseB[104] = 1250000.0*dwdx19;
    JSparseB[105] = -1250000.0*dwdx71;
    JSparseB[106] = 1250000.0*dwdx73 + 1250000.0*dwdx74 + 1250000.0*dwdx75 - 1250000.0*dwdx77;
}