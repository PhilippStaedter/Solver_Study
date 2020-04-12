#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_bachmann2(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -2.5*dwdx1;
    JSparse[1] = -2.5*dwdx0;
    JSparse[2] = 2.5*dwdx0;
    JSparse[3] = 2.5*dwdx3;
    JSparse[4] = -2.5*dwdx2;
    JSparse[5] = -3.6363636363636362*dwdx4;
    JSparse[6] = 3.6363636363636362*dwdx4;
    JSparse[7] = -3.6363636363636362*dwdx5;
    JSparse[8] = 3.6363636363636362*dwdx5;
    JSparse[9] = -3.6363636363636362*dwdx6;
    JSparse[10] = 3.6363636363636362*dwdx6;
    JSparse[11] = -3.6363636363636362*dwdx7;
    JSparse[12] = 3.6363636363636362*dwdx7;
    JSparse[13] = 2.5*dwdx8;
    JSparse[14] = -3.6363636363636362*dwdx8;
    JSparse[15] = -2.5*dwdx9;
    JSparse[16] = 2.5*dwdx9;
    JSparse[17] = -2.5*dwdx10;
    JSparse[18] = 2.5*dwdx10;
    JSparse[19] = -2.5*dwdx11;
    JSparse[20] = -2.5*dwdx12;
    JSparse[21] = 2.5*dwdx13;
    JSparse[22] = -2.5*dwdx13;
    JSparse[23] = 2.5*dwdx12;
    JSparse[24] = 2.5*dwdx16;
    JSparse[25] = -2.5*dwdx16 - 2.5*dwdx17 - 2.5*dwdx18;
    JSparse[26] = -2.5*dwdx14;
    JSparse[27] = 2.5*dwdx14;
    JSparse[28] = -2.5*dwdx15;
    JSparse[29] = 2.5*dwdx17;
    JSparse[30] = 2.5*dwdx18;
    JSparse[31] = 2.5*dwdx15;
    JSparse[32] = -2.5*dwdx19;
    JSparse[33] = 2.5*dwdx19;
    JSparse[34] = 2.5*dwdx21 + 2.5*dwdx22 + 2.5*dwdx23 + 2.5*dwdx24;
    JSparse[35] = -2.5*dwdx21;
    JSparse[36] = 2.5*dwdx20;
    JSparse[37] = -2.5*dwdx20;
    JSparse[38] = -2.5*dwdx24;
    JSparse[39] = -2.5*dwdx22;
    JSparse[40] = -2.5*dwdx23;
    JSparse[41] = -2.5*dwdx25;
    JSparse[42] = 2.5*dwdx25 - 2.5*dwdx28 - 2.5*dwdx30;
    JSparse[43] = -2.5*dwdx29;
    JSparse[44] = -2.5*dwdx26 - 2.5*dwdx27;
    JSparse[45] = 2.5*dwdx31 + 2.5*dwdx32;
    JSparse[46] = 2.5*dwdx28 - 2.5*dwdx31;
    JSparse[47] = 2.5*dwdx30 - 2.5*dwdx32;
    JSparse[48] = 2.5*dwdx26 + 2.5*dwdx27;
    JSparse[49] = 2.5*dwdx34;
    JSparse[50] = -2.5*dwdx33;
    JSparse[51] = -3.6363636363636362*dwdx35;
    JSparse[52] = 3.6363636363636362*dwdx35;
    JSparse[53] = -3.6363636363636362*dwdx36;
    JSparse[54] = 3.6363636363636362*dwdx36;
    JSparse[55] = -3.6363636363636362*dwdx37;
    JSparse[56] = 3.6363636363636362*dwdx37;
    JSparse[57] = -3.6363636363636362*dwdx38;
    JSparse[58] = 3.6363636363636362*dwdx38;
    JSparse[59] = 2.5*dwdx39;
    JSparse[60] = -3.6363636363636362*dwdx39;
    JSparse[61] = -2.5*dwdx40 - 2.5*dwdx41;
    JSparse[62] = 2.5*dwdx40 + 2.5*dwdx41;
    JSparse[63] = 3.6363636363636362*dwdx43;
    JSparse[64] = 3.6363636363636362*dwdx44;
    JSparse[65] = 2.5*dwdx42;
    JSparse[66] = -3.6363636363636362*dwdx42;
    JSparse[67] = 2.5*dwdx49;
    JSparse[68] = -2.5*dwdx45;
    JSparse[69] = -2.5*dwdx46;
    JSparse[70] = 2.5*dwdx46;
    JSparse[71] = -2.5*dwdx47 - 2.5*dwdx48;
    JSparse[72] = -2.5*dwdx49;
    JSparse[73] = 2.5*dwdx47 + 2.5*dwdx48;
    JSparse[74] = 2.5*dwdx55;
    JSparse[75] = -2.5*dwdx50;
    JSparse[76] = -2.5*dwdx51;
    JSparse[77] = 2.5*dwdx51;
    JSparse[78] = -2.5*dwdx52 - 2.5*dwdx53;
    JSparse[79] = 2.5*dwdx54;
    JSparse[80] = -2.5*dwdx54 - 2.5*dwdx55;
    JSparse[81] = 2.5*dwdx52 + 2.5*dwdx53;
    JSparse[82] = 2.5*dwdx59;
    JSparse[83] = -2.5*dwdx56;
    JSparse[84] = 2.5*dwdx56;
    JSparse[85] = -2.5*dwdx57;
    JSparse[86] = 2.5*dwdx58;
    JSparse[87] = -2.5*dwdx58 - 2.5*dwdx59;
    JSparse[88] = 2.5*dwdx57;
    JSparse[89] = 3.6363636363636362*dwdx60;
    JSparse[90] = -2.5*dwdx60;
}