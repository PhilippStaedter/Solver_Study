#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model29_beuke30(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -76923.076923076922*dwdx0 - 76923.076923076922*dwdx1 - 76923.076923076922*dwdx2 - 76923.076923076922*dwdx3;
    JSparse[1] = -76923.076923076922*dwdx0;
    JSparse[2] = 1250000.0*dwdx2;
    JSparse[3] = 76923.076923076922*dwdx0;
    JSparse[4] = 76923.076923076922*dwdx3;
    JSparse[5] = 76923.076923076922*dwdx5;
    JSparse[6] = -76923.076923076922*dwdx4;
    JSparse[7] = 76923.076923076922*dwdx7;
    JSparse[8] = -1250000.0*dwdx6 - 1250000.0*dwdx7;
    JSparse[9] = -1250000.0*dwdx8;
    JSparse[10] = 1250000.0*dwdx8;
    JSparse[11] = 1250000.0*dwdx9;
    JSparse[12] = -1250000.0*dwdx9;
    JSparse[13] = -76923.076923076922*dwdx10 - 76923.076923076922*dwdx11;
    JSparse[14] = 76923.076923076922*dwdx10 + 76923.076923076922*dwdx11;
    JSparse[15] = 76923.076923076922*dwdx12;
    JSparse[16] = -76923.076923076922*dwdx12;
    JSparse[17] = -386100.38610038609*dwdx13 - 386100.38610038609*dwdx14;
    JSparse[18] = 4255319.1489361702*dwdx15;
    JSparse[19] = 5714285.7142857146*dwdx16;
    JSparse[20] = -1250000.0*dwdx17;
    JSparse[21] = 1250000.0*dwdx17;
    JSparse[22] = 1250000.0*dwdx18;
    JSparse[23] = -1250000.0*dwdx18;
    JSparse[24] = 1250000.0*dwdx19 + 1250000.0*dwdx20;
    JSparse[25] = -1250000.0*dwdx20;
    JSparse[26] = -1250000.0*dwdx19;
    JSparse[27] = 386100.38610038609*dwdx21 - 386100.38610038609*dwdx22;
    JSparse[28] = -76923.076923076922*dwdx21;
    JSparse[29] = 386100.38610038609*dwdx22;
    JSparse[30] = -386100.38610038609*dwdx22;
    JSparse[31] = 386100.38610038609*dwdx23;
    JSparse[32] = -76923.076923076922*dwdx23;
    JSparse[33] = -76923.076923076922*dwdx25;
    JSparse[34] = 76923.076923076922*dwdx25;
    JSparse[35] = -386100.38610038609*dwdx26;
    JSparse[36] = 386100.38610038609*dwdx26 - 386100.38610038609*dwdx27;
    JSparse[37] = 76923.076923076922*dwdx27;
    JSparse[38] = -76923.076923076922*dwdx28;
    JSparse[39] = 76923.076923076922*dwdx28;
    JSparse[40] = -76923.076923076922*dwdx24;
    JSparse[41] = 76923.076923076922*dwdx24;
    JSparse[42] = -386100.38610038609*dwdx26;
    JSparse[43] = 76923.076923076922*dwdx29;
    JSparse[44] = -76923.076923076922*dwdx29;
    JSparse[45] = -2127659.5744680851*dwdx31;
    JSparse[46] = 386100.38610038609*dwdx30;
    JSparse[47] = -2857142.8571428573*dwdx33;
    JSparse[48] = 386100.38610038609*dwdx32;
    JSparse[49] = 2127659.5744680851*dwdx34;
    JSparse[50] = -2127659.5744680851*dwdx34;
    JSparse[51] = 2857142.8571428573*dwdx35;
    JSparse[52] = -2857142.8571428573*dwdx35;
    JSparse[53] = -76923.076923076922*dwdx36 - 76923.076923076922*dwdx37;
    JSparse[54] = 76923.076923076922*dwdx36 + 76923.076923076922*dwdx37;
    JSparse[55] = -1250000.0*dwdx38;
    JSparse[56] = 1250000.0*dwdx38;
    JSparse[57] = 76923.076923076922*dwdx39;
    JSparse[58] = -76923.076923076922*dwdx39;
    JSparse[59] = 1250000.0*dwdx42;
    JSparse[60] = -1250000.0*dwdx40 - 1250000.0*dwdx41;
    JSparse[61] = 1250000.0*dwdx41;
    JSparse[62] = -1250000.0*dwdx43 - 1250000.0*dwdx44 - 1250000.0*dwdx45 + 1250000.0*dwdx46;
    JSparse[63] = 1250000.0*dwdx43 - 1250000.0*dwdx46;
    JSparse[64] = 1250000.0*dwdx44 - 1250000.0*dwdx46;
    JSparse[65] = 76923.076923076922*dwdx45;
    JSparse[66] = -76923.076923076922*dwdx47;
    JSparse[67] = -76923.076923076922*dwdx47 - 76923.076923076922*dwdx48;
    JSparse[68] = 1250000.0*dwdx48;
    JSparse[69] = 76923.076923076922*dwdx47;
    JSparse[70] = 76923.076923076922*dwdx49;
    JSparse[71] = 1250000.0*dwdx52;
    JSparse[72] = 1250000.0*dwdx51;
    JSparse[73] = 1250000.0*dwdx53;
    JSparse[74] = -76923.076923076922*dwdx50;
    JSparse[75] = 1250000.0*dwdx50 - 1250000.0*dwdx51 - 1250000.0*dwdx53;
    JSparse[76] = -1250000.0*dwdx53;
    JSparse[77] = -76923.076923076922*dwdx55;
    JSparse[78] = 1250000.0*dwdx56;
    JSparse[79] = -1250000.0*dwdx56;
    JSparse[80] = -1250000.0*dwdx54 + 1250000.0*dwdx55 - 1250000.0*dwdx56;
    JSparse[81] = -76923.076923076922*dwdx57 - 76923.076923076922*dwdx58;
    JSparse[82] = 76923.076923076922*dwdx57 + 76923.076923076922*dwdx58;
    JSparse[83] = -76923.076923076922*dwdx59;
    JSparse[84] = 76923.076923076922*dwdx61;
    JSparse[85] = -76923.076923076922*dwdx61;
    JSparse[86] = -76923.076923076922*dwdx60;
    JSparse[87] = 76923.076923076922*dwdx60;
    JSparse[88] = 76923.076923076922*dwdx59 + 76923.076923076922*dwdx60;
    JSparse[89] = -386100.38610038609*dwdx62;
    JSparse[90] = 386100.38610038609*dwdx62;
    JSparse[91] = -386100.38610038609*dwdx62 - 386100.38610038609*dwdx63;
    JSparse[92] = -76923.076923076922*dwdx64 + 76923.076923076922*dwdx66;
    JSparse[93] = -1250000.0*dwdx67;
    JSparse[94] = -76923.076923076922*dwdx64 + 76923.076923076922*dwdx65;
    JSparse[95] = 76923.076923076922*dwdx64 - 76923.076923076922*dwdx65 - 76923.076923076922*dwdx66 + 76923.076923076922*dwdx67 - 76923.076923076922*dwdx68;
    JSparse[96] = 76923.076923076922*dwdx68;
    JSparse[97] = 76923.076923076922*dwdx68;
    JSparse[98] = 76923.076923076922*dwdx70;
    JSparse[99] = -76923.076923076922*dwdx69 - 76923.076923076922*dwdx70 - 76923.076923076922*dwdx71;
    JSparse[100] = 1250000.0*dwdx71;
    JSparse[101] = -76923.076923076922*dwdx72;
    JSparse[102] = 1250000.0*dwdx76;
    JSparse[103] = 1250000.0*dwdx75;
    JSparse[104] = 1250000.0*dwdx74;
    JSparse[105] = -76923.076923076922*dwdx77;
    JSparse[106] = -1250000.0*dwdx73 - 1250000.0*dwdx74 - 1250000.0*dwdx75 + 1250000.0*dwdx77;
}