#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_schilling1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0;
    JSparse[1] = 1.0*dwdx0;
    JSparse[2] = -1.0*dwdx1;
    JSparse[3] = 1.0*dwdx1;
    JSparse[4] = -1.0*dwdx2;
    JSparse[5] = 1.0*dwdx2;
    JSparse[6] = -1.0*dwdx3;
    JSparse[7] = 1.0*dwdx3;
    JSparse[8] = -1.0*dwdx4;
    JSparse[9] = 1.0*dwdx4;
    JSparse[10] = -1.0*dwdx5;
    JSparse[11] = 1.0*dwdx5;
    JSparse[12] = -1.0*dwdx6;
    JSparse[13] = 1.0*dwdx6;
    JSparse[14] = -1.0*dwdx7;
    JSparse[15] = 1.0*dwdx7;
    JSparse[16] = -1.0*dwdx8 - 1.0*dwdx9;
    JSparse[17] = 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[18] = -1.0*dwdx10 - 1.0*dwdx11;
    JSparse[19] = 1.0*dwdx10 + 1.0*dwdx11;
    JSparse[20] = -1.0*dwdx12;
    JSparse[21] = 1.0*dwdx12;
    JSparse[22] = -1.0*dwdx13;
    JSparse[23] = 1.0*dwdx13;
    JSparse[24] = -1.0*dwdx14;
    JSparse[25] = 1.0*dwdx14;
    JSparse[26] = -1.0*dwdx15;
    JSparse[27] = 1.0*dwdx15;
    JSparse[28] = -1.0*dwdx16;
    JSparse[29] = 1.0*dwdx16;
    JSparse[30] = -1.0*dwdx17;
    JSparse[31] = 1.0*dwdx17;
    JSparse[32] = -1.0*dwdx18;
    JSparse[33] = 1.0*dwdx18;
    JSparse[34] = -1.0*dwdx19;
    JSparse[35] = 1.0*dwdx19;
    JSparse[36] = 1.0*dwdx21;
    JSparse[37] = 1.0*dwdx22;
    JSparse[38] = 1.0*dwdx20;
    JSparse[39] = -1.0*dwdx20;
    JSparse[40] = -1.0*dwdx21;
    JSparse[41] = -1.0*dwdx22;
    JSparse[42] = 1.0*dwdx23;
    JSparse[43] = -1.0*dwdx23;
    JSparse[44] = -1.0*dwdx25;
    JSparse[45] = 1.0*dwdx24;
    JSparse[46] = -1.0*dwdx24 - 1.0*dwdx26 - 1.0*dwdx27;
    JSparse[47] = 1.0*dwdx25;
    JSparse[48] = 1.0*dwdx26 + 1.0*dwdx27;
    JSparse[49] = 1.0*dwdx30;
    JSparse[50] = -1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30;
    JSparse[51] = 1.0*dwdx28 + 1.0*dwdx29;
    JSparse[52] = 1.0*dwdx33;
    JSparse[53] = -1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33;
    JSparse[54] = 1.0*dwdx31 + 1.0*dwdx32;
    JSparse[55] = 1.0*dwdx34;
    JSparse[56] = -1.0*dwdx36;
    JSparse[57] = -1.0*dwdx35;
    JSparse[58] = 1.0*dwdx36;
    JSparse[59] = 1.0*dwdx35;
    JSparse[60] = -1.0*dwdx34;
    JSparse[61] = -1.0*dwdx38;
    JSparse[62] = 1.0*dwdx37;
    JSparse[63] = 1.0*dwdx38;
    JSparse[64] = -1.0*dwdx37;
    JSparse[65] = 1.0*dwdx40;
    JSparse[66] = -1.0*dwdx39 - 1.0*dwdx40;
    JSparse[67] = 1.0*dwdx39;
    JSparse[68] = 1.0*dwdx42;
    JSparse[69] = -1.0*dwdx41 - 1.0*dwdx42;
    JSparse[70] = 1.0*dwdx41;
    JSparse[71] = -1.0*dwdx45;
    JSparse[72] = -1.0*dwdx44;
    JSparse[73] = 1.0*dwdx43;
    JSparse[74] = 1.0*dwdx45 - 1.0*dwdx47;
    JSparse[75] = 1.0*dwdx44 - 1.0*dwdx46;
    JSparse[76] = -1.0*dwdx43;
    JSparse[77] = 1.0*dwdx47;
    JSparse[78] = 1.0*dwdx46;
    JSparse[79] = 1.0*dwdx48;
    JSparse[80] = -1.0*dwdx48;
    JSparse[81] = -1.0*dwdx50;
    JSparse[82] = 1.0*dwdx49;
    JSparse[83] = 1.0*dwdx50;
    JSparse[84] = -1.0*dwdx49;
    JSparse[85] = -1.0*dwdx52;
    JSparse[86] = 1.0*dwdx51;
    JSparse[87] = 1.0*dwdx52;
    JSparse[88] = -1.0*dwdx51;
    JSparse[89] = -1.0*dwdx54;
    JSparse[90] = -1.0*dwdx55;
    JSparse[91] = 1.0*dwdx54 - 1.0*dwdx56;
    JSparse[92] = 1.0*dwdx55 - 1.0*dwdx57;
    JSparse[93] = 1.0*dwdx53;
    JSparse[94] = 1.0*dwdx56;
    JSparse[95] = 1.0*dwdx57;
    JSparse[96] = -1.0*dwdx53;
    JSparse[97] = -1.0*dwdx59;
    JSparse[98] = -1.0*dwdx60;
    JSparse[99] = 1.0*dwdx59 - 1.0*dwdx61;
    JSparse[100] = 1.0*dwdx60 - 1.0*dwdx62;
    JSparse[101] = 1.0*dwdx58;
    JSparse[102] = 1.0*dwdx61;
    JSparse[103] = 1.0*dwdx62;
    JSparse[104] = -1.0*dwdx58;
}