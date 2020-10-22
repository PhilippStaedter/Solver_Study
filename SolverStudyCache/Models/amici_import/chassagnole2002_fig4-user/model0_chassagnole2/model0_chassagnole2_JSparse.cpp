#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_chassagnole2(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    JSparse[1] = -1.0*dwdx0;
    JSparse[2] = 1.0*dwdx0 + 1.0*dwdx3;
    JSparse[3] = -1.0*dwdx4 - 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7;
    JSparse[4] = 1.0*dwdx6 + 1.0*dwdx7;
    JSparse[5] = -1.0*dwdx6 + 1.0*dwdx7;
    JSparse[6] = -1.0*dwdx4;
    JSparse[7] = -1.0*dwdx6;
    JSparse[8] = -1.0*dwdx7;
    JSparse[9] = 1.0*dwdx10 - 1.0*dwdx11;
    JSparse[10] = 1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[11] = 1.0*dwdx8;
    JSparse[12] = -1.0*dwdx9;
    JSparse[13] = -1.0*dwdx10 + 1.0*dwdx11;
    JSparse[14] = -1.0*dwdx10;
    JSparse[15] = -1.0*dwdx11;
    JSparse[16] = 1.0*dwdx13;
    JSparse[17] = -1.0*dwdx13 - 1.0*dwdx16;
    JSparse[18] = -1.0*dwdx14;
    JSparse[19] = 1.0*dwdx13;
    JSparse[20] = -1.0*dwdx15 - 1.0*dwdx17;
    JSparse[21] = 1.0*dwdx15;
    JSparse[22] = -1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
    JSparse[23] = -1.0*dwdx20;
    JSparse[24] = 1.0*dwdx23;
    JSparse[25] = 1.0*dwdx24;
    JSparse[26] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 + 65.0*dwdx25;
    JSparse[27] = -1.0*dwdx25;
    JSparse[28] = -65.0*dwdx25;
    JSparse[29] = 1.0*dwdx22;
    JSparse[30] = 65.0*dwdx25;
    JSparse[31] = 1.0*dwdx26 - 1.0*dwdx30;
    JSparse[32] = 1.0*dwdx29 - 1.0*dwdx32;
    JSparse[33] = 1.0*dwdx29 + 1.0*dwdx32;
    JSparse[34] = -1.0*dwdx26;
    JSparse[35] = 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29 + 1.0*dwdx30 + 1.0*dwdx31 + 1.0*dwdx32;
    JSparse[36] = 1.0*dwdx28;
    JSparse[37] = -1.0*dwdx31;
    JSparse[38] = -1.0*dwdx29 + 1.0*dwdx31;
    JSparse[39] = -1.0*dwdx31 - 1.0*dwdx32;
    JSparse[40] = 65.0*dwdx34;
    JSparse[41] = 1.0*dwdx33 - 1.0*dwdx34;
    JSparse[42] = -65.0*dwdx34;
    JSparse[43] = 65.0*dwdx34;
    JSparse[44] = -1.0*dwdx35;
    JSparse[45] = -1.0*dwdx38;
    JSparse[46] = 1.0*dwdx38;
    JSparse[47] = 65.0*dwdx40;
    JSparse[48] = -1.0*dwdx40;
    JSparse[49] = -1.0*dwdx35 + 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx39 - 65.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42;
    JSparse[50] = -1.0*dwdx36;
    JSparse[51] = 1.0*dwdx39 + 65.0*dwdx40;
    JSparse[52] = 1.0*dwdx45;
    JSparse[53] = -1.0*dwdx45;
    JSparse[54] = -1.0*dwdx43 - 1.0*dwdx44;
    JSparse[55] = 1.0*dwdx44;
    JSparse[56] = 1.0*dwdx46;
    JSparse[57] = -1.0*dwdx46 - 1.0*dwdx47 + 1.0*dwdx48;
    JSparse[58] = -1.0*dwdx48;
    JSparse[59] = 1.0*dwdx51;
    JSparse[60] = -1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    JSparse[61] = -1.0*dwdx50;
    JSparse[62] = -1.0*dwdx53;
    JSparse[63] = 1.0*dwdx54;
    JSparse[64] = 1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55;
    JSparse[65] = 65.0*dwdx57;
    JSparse[66] = -1.0*dwdx57;
    JSparse[67] = -65.0*dwdx57;
    JSparse[68] = -1.0*dwdx56 + 65.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    JSparse[69] = 1.0*dwdx63;
    JSparse[70] = -1.0*dwdx60 + 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63;
    JSparse[71] = -1.0*dwdx61;
    JSparse[72] = 1.0*dwdx63;
    JSparse[73] = -1.0*dwdx63;
    JSparse[74] = 1.0*dwdx64;
    JSparse[75] = -1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66;
    JSparse[76] = 1.0*dwdx66;
    JSparse[77] = 1.0*dwdx68;
    JSparse[78] = 1.0*dwdx68;
    JSparse[79] = -1.0*dwdx68 + 1.0*dwdx69;
    JSparse[80] = -1.0*dwdx69;
    JSparse[81] = -1.0*dwdx67 - 1.0*dwdx68 + 1.0*dwdx69;
    JSparse[82] = -1.0*dwdx69;
    JSparse[83] = -1.0*dwdx72;
    JSparse[84] = 1.0*dwdx72;
    JSparse[85] = 1.0*dwdx71 + 1.0*dwdx72;
    JSparse[86] = -1.0*dwdx71;
    JSparse[87] = -1.0*dwdx70;
    JSparse[88] = 1.0*dwdx71;
    JSparse[89] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73;
}