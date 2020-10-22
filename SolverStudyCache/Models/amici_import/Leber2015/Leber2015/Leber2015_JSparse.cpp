#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_Leber2015(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx3 - 1.0*dwdx7;
    JSparse[1] = 1.0*dwdx7;
    JSparse[2] = 1.0*dwdx2;
    JSparse[3] = -0.25*dwdx5;
    JSparse[4] = 0.25*dwdx5;
    JSparse[5] = 0.25*dwdx6;
    JSparse[6] = 14.285714285714285*dwdx0;
    JSparse[7] = -14.285714285714285*dwdx4;
    JSparse[8] = 14.285714285714285*dwdx4;
    JSparse[9] = -1.0*dwdx10;
    JSparse[10] = -1.0*dwdx9;
    JSparse[11] = 1.0*dwdx9;
    JSparse[12] = 1.0*dwdx10;
    JSparse[13] = -1.0*dwdx8;
    JSparse[14] = -1.0*dwdx11;
    JSparse[15] = 1.0*dwdx11;
    JSparse[16] = -1.0*dwdx13;
    JSparse[17] = -1.0*dwdx12;
    JSparse[18] = 1.0*dwdx12 - 1.0*dwdx15;
    JSparse[19] = 1.0*dwdx13;
    JSparse[20] = -1.0*dwdx14;
    JSparse[21] = 1.0*dwdx14;
    JSparse[22] = -1.0*dwdx16;
    JSparse[23] = 1.0*dwdx16;
    JSparse[24] = -1.0*dwdx17 - 1.0*dwdx18;
    JSparse[25] = 1.0*dwdx18;
    JSparse[26] = -1.0*dwdx19 + 1.0*dwdx20;
    JSparse[27] = -1.0*dwdx21;
    JSparse[28] = -1.0*dwdx24;
    JSparse[29] = -1.0*dwdx26;
    JSparse[30] = 1.0*dwdx26;
    JSparse[31] = -1.0*dwdx22;
    JSparse[32] = -0.25*dwdx23;
    JSparse[33] = 0.25*dwdx23 + 0.25*dwdx25;
    JSparse[34] = -0.25*dwdx25;
    JSparse[35] = -1.0*dwdx29;
    JSparse[36] = 1.0*dwdx29;
    JSparse[37] = -0.25*dwdx27 - 0.25*dwdx28;
    JSparse[38] = 0.25*dwdx27;
    JSparse[39] = 0.25*dwdx28;
    JSparse[40] = -1.0*dwdx30;
    JSparse[41] = 1.0*dwdx30;
    JSparse[42] = 1.0*dwdx31;
    JSparse[43] = 0.25*dwdx32;
    JSparse[44] = -0.25*dwdx32;
    JSparse[45] = -1.0*dwdx35;
    JSparse[46] = -1.0*dwdx34;
    JSparse[47] = 1.0*dwdx34;
    JSparse[48] = 1.0*dwdx35;
    JSparse[49] = -1.0*dwdx36;
    JSparse[50] = 0.25*dwdx33 + 0.25*dwdx37;
    JSparse[51] = -0.25*dwdx33 - 0.25*dwdx37;
    JSparse[52] = -1.0*dwdx39;
    JSparse[53] = -0.25*dwdx38;
    JSparse[54] = 0.25*dwdx38 + 0.25*dwdx40;
    JSparse[55] = -0.25*dwdx40;
    JSparse[56] = -0.25*dwdx41;
    JSparse[57] = -14.285714285714285*dwdx42;
    JSparse[58] = 1.0*dwdx42;
    JSparse[59] = 0.25*dwdx43;
    JSparse[60] = -1.0*dwdx45;
    JSparse[61] = 1.0*dwdx44;
    JSparse[62] = 1.0*dwdx48;
    JSparse[63] = -0.25*dwdx47;
    JSparse[64] = 0.25*dwdx47 + 0.25*dwdx50;
    JSparse[65] = -0.25*dwdx50;
    JSparse[66] = 0.25*dwdx51;
    JSparse[67] = -14.285714285714285*dwdx46 - 14.285714285714285*dwdx49;
    JSparse[68] = 14.285714285714285*dwdx49;
    JSparse[69] = -14.285714285714285*dwdx52;
    JSparse[70] = 1.0*dwdx54;
    JSparse[71] = 0.25*dwdx56;
    JSparse[72] = -14.285714285714285*dwdx55;
    JSparse[73] = -14.285714285714285*dwdx53 + 14.285714285714285*dwdx55;
    JSparse[74] = -1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    JSparse[75] = 1.0*dwdx58;
    JSparse[76] = 1.0*dwdx59;
    JSparse[77] = 14.285714285714285*dwdx60;
    JSparse[78] = -1.0*dwdx60;
    JSparse[79] = 14.285714285714285*dwdx61;
    JSparse[80] = -1.0*dwdx61;
    JSparse[81] = 14.285714285714285*dwdx62;
    JSparse[82] = -1.0*dwdx62;
}